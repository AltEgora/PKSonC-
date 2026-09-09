#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace {
    using ui8 = uint8_t;
    using ui32 = uint32_t;
    using ui64 = uint64_t;
    using i64 = int64_t;
    using ui128 = unsigned __int128;

    const int NUMERATORS[] = {1, 1, 1, 1, 1, 3, 1, 5, 3, 7, 1, 5, 3, 7, 1, 3, 2, 3, 4, 8};
    const int DENOMINATORS[] = {32, 16, 8, 6, 4, 8, 2, 8, 4, 8, 1, 4, 2, 2, 1, 1, 1, 1, 1, 1};

    class TGraphCodecError: public std::runtime_error {
    public:
        explicit TGraphCodecError(const std::string& message)
            : std::runtime_error(message)
        {
        }
    };

    std::vector<char> ReadWholeFile(const std::string& path) {
        FILE* file = std::fopen(path.c_str(), "rb");
        if (file == nullptr) {
            throw TGraphCodecError("cannot open input file: " + path);
        }
        std::fseek(file, 0, SEEK_END);
        const long size = std::ftell(file);
        if (size < 0) {
            std::fclose(file);
            throw TGraphCodecError("cannot determine size of: " + path);
        }
        std::fseek(file, 0, SEEK_SET);
        std::vector<char> buffer(static_cast<size_t>(size));
        const size_t read = std::fread(buffer.data(), 1, static_cast<size_t>(size), file);
        std::fclose(file);
        if (read != static_cast<size_t>(size)) {
            throw TGraphCodecError("failed to read: " + path);
        }
        return buffer;
    }

    inline int FloorLog2(ui64 value) {
        return 63 - __builtin_clzll(value);
    }

    ui64 RowStart(ui64 row, ui64 vertexCount) {
        if (row == 0) {
            return 0;
        }
        const ui128 rank = static_cast<ui128>(row) * vertexCount - static_cast<ui128>(row) * (row - 1) / 2;
        return static_cast<ui64>(rank);
    }

    class TBitWriter {
    public:
        void Write(ui64 value, int bitCount) {
            if (bitCount == 0) {
                return;
            }
            Acc_ = (Acc_ << bitCount) | (value & ((1ULL << bitCount) - 1));
            AccBits_ += bitCount;
            BitLength_ += static_cast<ui64>(bitCount);
            while (AccBits_ >= 8) {
                AccBits_ -= 8;
                Bytes_.push_back(static_cast<ui8>((Acc_ >> AccBits_) & 0xFF));
            }
        }

        void WriteUnary(ui64 count) {
            while (count >= 32) {
                Write(~0ULL >> 32, 32);
                count -= 32;
            }
            Write((1ULL << count) - 1, static_cast<int>(count));
            Write(0, 1);
        }

        void Finish() {
            if (AccBits_ > 0) {
                Bytes_.push_back(static_cast<ui8>((Acc_ << (8 - AccBits_)) & 0xFF));
                AccBits_ = 0;
            }
        }

        ui64 BitLength() const {
            return BitLength_;
        }

        const std::vector<ui8>& Bytes() const {
            return Bytes_;
        }

    private:
        std::vector<ui8> Bytes_;
        ui64 Acc_ = 0;
        int AccBits_ = 0;
        ui64 BitLength_ = 0;
    };

    class TBitReader {
    public:
        TBitReader(const ui8* data, size_t size)
            : Data_(data)
            , Size_(size)
        {
        }

        ui64 Read(int bitCount) {
            if (bitCount == 0) {
                return 0;
            }
            while (AccBits_ < bitCount) {
                const ui8 byte = (Pos_ < Size_) ? Data_[Pos_] : 0;
                ++Pos_;
                Acc_ = (Acc_ << 8) | byte;
                AccBits_ += 8;
            }
            AccBits_ -= bitCount;
            return (Acc_ >> AccBits_) & ((1ULL << bitCount) - 1);
        }

        ui64 ReadUnary() {
            ui64 count = 0;
            while (Read(1) == 1) {
                ++count;
            }
            return count;
        }

    private:
        const ui8* Data_;
        size_t Size_;
        size_t Pos_ = 0;
        ui64 Acc_ = 0;
        int AccBits_ = 0;
    };

    inline void WriteMinBinary(TBitWriter* writer, ui64 value, ui64 range) {
        if (range <= 1) {
            return;
        }
        const int bits = FloorLog2(range);
        const ui64 shortCount = (1ULL << (bits + 1)) - range;
        if (value < shortCount) {
            writer->Write(value, bits);
        } else {
            writer->Write(value + shortCount, bits + 1);
        }
    }

    inline ui64 ReadMinBinary(TBitReader* reader, ui64 range) {
        if (range <= 1) {
            return 0;
        }
        const int bits = FloorLog2(range);
        const ui64 shortCount = (1ULL << (bits + 1)) - range;
        ui64 value = reader->Read(bits);
        if (value >= shortCount) {
            value = ((value << 1) | reader->Read(1)) - shortCount;
        }
        return value;
    }

    inline ui64 CenterBase(ui64 range) {
        const int bits = FloorLog2(range);
        const ui64 shortCount = (1ULL << (bits + 1)) - range;
        return (range - shortCount) / 2;
    }

    inline void WriteCenteredMinBinary(TBitWriter* writer, ui64 value, ui64 range) {
        if (range <= 1) {
            return;
        }
        const ui64 base = CenterBase(range);
        const ui64 rotated = (value >= base) ? (value - base) : (value + range - base);
        WriteMinBinary(writer, rotated, range);
    }

    inline ui64 ReadCenteredMinBinary(TBitReader* reader, ui64 range) {
        if (range <= 1) {
            return 0;
        }
        ui64 value = ReadMinBinary(reader, range) + CenterBase(range);
        if (value >= range) {
            value -= range;
        }
        return value;
    }

    ui64 GolombCost(const std::vector<ui64>& gaps, ui64 parameter) {
        const int bits = FloorLog2(parameter);
        const ui64 shortCount = (1ULL << (bits + 1)) - parameter;
        ui64 total = 0;
        for (const ui64 gap : gaps) {
            const ui64 quotient = gap / parameter;
            const ui64 remainder = gap - quotient * parameter;
            total += quotient + 1;
            if (parameter > 1) {
                total += (remainder < shortCount) ? static_cast<ui64>(bits) : static_cast<ui64>(bits + 1);
            }
        }
        return total;
    }

    ui64 ChooseGolombParameter(const std::vector<ui64>& gaps) {
        if (gaps.empty()) {
            return 1;
        }

        ui64 sum = 0;
        ui64 maxGap = 0;
        for (const ui64 gap : gaps) {
            sum += gap;
            maxGap = std::max(maxGap, gap);
        }
        const ui64 mean = sum / gaps.size() + 1;

        std::vector<ui64> candidates;
        candidates.reserve(sizeof(NUMERATORS) / sizeof(NUMERATORS[0]) + 2);
        for (size_t index = 0; index < sizeof(NUMERATORS) / sizeof(NUMERATORS[0]); ++index) {
            const ui64 candidate = mean * static_cast<ui64>(NUMERATORS[index]) / static_cast<ui64>(DENOMINATORS[index]);
            candidates.push_back(std::max<ui64>(candidate, 1));
        }
        candidates.push_back(1);
        candidates.push_back(std::max<ui64>(maxGap, 1));

        ui64 bestParameter = 1;
        ui64 bestBits = ~0ULL;
        for (const ui64 candidate : candidates) {
            const ui64 bits = GolombCost(gaps, candidate);
            if (bits < bestBits) {
                bestBits = bits;
                bestParameter = candidate;
            }
        }
        return bestParameter;
    }

    void WriteGolomb(TBitWriter* writer, ui64 value, ui64 parameter) {
        const ui64 quotient = value / parameter;
        writer->WriteUnary(quotient);
        WriteMinBinary(writer, value - quotient * parameter, parameter);
    }

    ui64 ReadGolomb(TBitReader* reader, ui64 parameter) {
        const ui64 quotient = reader->ReadUnary();
        return quotient * parameter + ReadMinBinary(reader, parameter);
    }

    void BicEncode(const ui64* values, i64 left, i64 right, ui64 low, ui64 high, TBitWriter* writer) {
        if (left > right) {
            return;
        }
        const i64 middle = left + (right - left) / 2;
        const ui64 rangeLow = low + static_cast<ui64>(middle - left);
        const ui64 rangeHigh = high - static_cast<ui64>(right - middle);
        const ui64 value = values[middle];
        WriteCenteredMinBinary(writer, value - rangeLow, rangeHigh - rangeLow + 1);
        if (middle > left) {
            BicEncode(values, left, middle - 1, low, value - 1, writer);
        }
        if (middle < right) {
            BicEncode(values, middle + 1, right, value + 1, high, writer);
        }
    }

    void BicDecode(ui64* values, i64 left, i64 right, ui64 low, ui64 high, TBitReader* reader) {
        if (left > right) {
            return;
        }
        const i64 middle = left + (right - left) / 2;
        const ui64 rangeLow = low + static_cast<ui64>(middle - left);
        const ui64 rangeHigh = high - static_cast<ui64>(right - middle);
        const ui64 value = rangeLow + ReadCenteredMinBinary(reader, rangeHigh - rangeLow + 1);
        values[middle] = value;
        if (middle > left) {
            BicDecode(values, left, middle - 1, low, value - 1, reader);
        }
        if (middle < right) {
            BicDecode(values, middle + 1, right, value + 1, high, reader);
        }
    }

    void WriteUi64(FILE* file, ui64 value) {
        ui8 buffer[8];
        for (int index = 0; index < 8; ++index) {
            buffer[index] = static_cast<ui8>(value & 0xFF);
            value >>= 8;
        }
        std::fwrite(buffer, 1, 8, file);
    }

    ui64 ReadUi64(const char* source) {
        ui64 value = 0;
        for (int index = 7; index >= 0; --index) {
            value = (value << 8) | static_cast<ui8>(source[index]);
        }
        return value;
    }

    inline const char* ParseUInt(const char* position, const char* end, ui32* result) {
        ui64 value = 0;
        while (position < end && *position >= '0' && *position <= '9') {
            value = value * 10 + static_cast<ui64>(*position - '0');
            ++position;
        }
        *result = static_cast<ui32>(value);
        return position;
    }

    void ParseTsv(const std::vector<char>& buffer, std::vector<ui32>* firstIds, std::vector<ui32>* secondIds,
                  std::vector<ui8>* weights) {
        const char* position = buffer.data();
        const char* end = position + buffer.size();
        while (position < end) {
            while (position < end && (*position == '\n' || *position == '\r')) {
                ++position;
            }
            if (position >= end) {
                break;
            }
            ui32 first = 0;
            ui32 second = 0;
            ui32 weight = 0;
            position = ParseUInt(position, end, &first);
            if (position < end && *position == '\t') {
                ++position;
            }
            position = ParseUInt(position, end, &second);
            if (position < end && *position == '\t') {
                ++position;
            }
            position = ParseUInt(position, end, &weight);
            firstIds->push_back(first);
            secondIds->push_back(second);
            weights->push_back(static_cast<ui8>(weight));
            while (position < end && *position != '\n') {
                ++position;
            }
            if (position < end) {
                ++position;
            }
        }
    }

    inline void AppendUInt(std::string* output, ui32 value) {
        if (value == 0) {
            output->push_back('0');
            return;
        }
        char digits[10];
        int count = 0;
        while (value > 0) {
            digits[count++] = static_cast<char>('0' + value % 10);
            value /= 10;
        }
        for (int index = count - 1; index >= 0; --index) {
            output->push_back(digits[index]);
        }
    }

    void Serialize(const std::string& inputPath, const std::string& outputPath) {
        std::vector<char> buffer = ReadWholeFile(inputPath);

        std::vector<ui32> firstIds;
        std::vector<ui32> secondIds;
        std::vector<ui8> weights;
        const size_t estimate = buffer.size() / 12 + 16;
        firstIds.reserve(estimate);
        secondIds.reserve(estimate);
        weights.reserve(estimate);
        ParseTsv(buffer, &firstIds, &secondIds, &weights);
        std::vector<char>().swap(buffer);

        const ui64 edgeCount = firstIds.size();

        std::vector<ui32> verts;
        verts.reserve(2 * edgeCount);
        verts.insert(verts.end(), firstIds.begin(), firstIds.end());
        verts.insert(verts.end(), secondIds.begin(), secondIds.end());
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        const ui64 vertexCount = verts.size();

        std::vector<ui64> packedEdges(edgeCount);
        for (ui64 index = 0; index < edgeCount; ++index) {
            const ui64 firstIndex =
                static_cast<ui64>(std::lower_bound(verts.begin(), verts.end(), firstIds[index]) - verts.begin());
            const ui64 secondIndex =
                static_cast<ui64>(std::lower_bound(verts.begin(), verts.end(), secondIds[index]) - verts.begin());
            const ui64 row = std::min(firstIndex, secondIndex);
            const ui64 column = std::max(firstIndex, secondIndex);
            packedEdges[index] = ((RowStart(row, vertexCount) + (column - row)) << 8) | weights[index];
        }
        std::vector<ui32>().swap(firstIds);
        std::vector<ui32>().swap(secondIds);
        std::vector<ui8>().swap(weights);

        std::sort(packedEdges.begin(), packedEdges.end());

        std::vector<ui64> ranks(edgeCount);
        std::vector<ui8> sortedWeights(edgeCount);
        for (ui64 index = 0; index < edgeCount; ++index) {
            ranks[index] = packedEdges[index] >> 8;
            sortedWeights[index] = static_cast<ui8>(packedEdges[index] & 0xFF);
        }
        std::vector<ui64>().swap(packedEdges);

        std::vector<ui64> gaps(vertexCount);
        for (ui64 index = 0; index < vertexCount; ++index) {
            gaps[index] = (index == 0) ? verts[index] : (verts[index] - verts[index - 1] - 1);
        }
        const ui64 golombM = ChooseGolombParameter(gaps);
        TBitWriter idWriter;
        for (const ui64 gap : gaps) {
            WriteGolomb(&idWriter, gap, golombM);
        }
        const ui64 idBits = idWriter.BitLength();
        idWriter.Finish();
        std::vector<ui64>().swap(gaps);
        std::vector<ui32>().swap(verts);

        TBitWriter edgeWriter;
        if (edgeCount > 0) {
            const ui64 pairUniverse = RowStart(vertexCount, vertexCount);
            BicEncode(ranks.data(), 0, static_cast<i64>(edgeCount) - 1, 0, pairUniverse - 1, &edgeWriter);
        }
        const ui64 edgeBits = edgeWriter.BitLength();
        edgeWriter.Finish();
        std::vector<ui64>().swap(ranks);

        FILE* output = std::fopen(outputPath.c_str(), "wb");
        if (output == nullptr) {
            throw TGraphCodecError("cannot open output file: " + outputPath);
        }
        std::fwrite("GCP2", 1, 4, output);
        WriteUi64(output, vertexCount);
        WriteUi64(output, edgeCount);
        WriteUi64(output, golombM);
        WriteUi64(output, idBits);
        WriteUi64(output, edgeBits);
        if (!idWriter.Bytes().empty()) {
            std::fwrite(idWriter.Bytes().data(), 1, idWriter.Bytes().size(), output);
        }
        if (!edgeWriter.Bytes().empty()) {
            std::fwrite(edgeWriter.Bytes().data(), 1, edgeWriter.Bytes().size(), output);
        }
        if (!sortedWeights.empty()) {
            std::fwrite(sortedWeights.data(), 1, sortedWeights.size(), output);
        }
        std::fclose(output);
    }

    void Deserialize(const std::string& inputPath, const std::string& outputPath) {
        std::vector<char> buffer = ReadWholeFile(inputPath);
        if (buffer.size() < 44 || std::memcmp(buffer.data(), "GCP2", 4) != 0) {
            throw TGraphCodecError("bad magic in: " + inputPath);
        }

        const char* position = buffer.data() + 4;
        const ui64 vertexCount = ReadUi64(position);
        const ui64 edgeCount = ReadUi64(position + 8);
        const ui64 golombM = ReadUi64(position + 16);
        const ui64 idBits = ReadUi64(position + 24);
        const ui64 edgeBits = ReadUi64(position + 32);
        position += 40;

        const ui64 idByteCount = (idBits + 7) / 8;
        const ui64 edgeByteCount = (edgeBits + 7) / 8;

        std::vector<ui32> verts(vertexCount);
        TBitReader idReader(reinterpret_cast<const ui8*>(position), idByteCount);
        ui64 previous = 0;
        for (ui64 index = 0; index < vertexCount; ++index) {
            const ui64 gap = ReadGolomb(&idReader, golombM);
            const ui64 id = (index == 0) ? gap : (previous + gap + 1);
            verts[index] = static_cast<ui32>(id);
            previous = id;
        }
        position += idByteCount;

        std::vector<ui64> ranks(edgeCount);
        if (edgeCount > 0) {
            TBitReader edgeReader(reinterpret_cast<const ui8*>(position), edgeByteCount);
            const ui64 pairUniverse = RowStart(vertexCount, vertexCount);
            BicDecode(ranks.data(), 0, static_cast<i64>(edgeCount) - 1, 0, pairUniverse - 1, &edgeReader);
        }
        position += edgeByteCount;

        const ui8* weights = reinterpret_cast<const ui8*>(position);

        std::string text;
        text.reserve(edgeCount * 24 + 16);

        ui64 row = 0;
        ui64 rowStart = 0;
        ui64 rowLength = vertexCount;
        for (ui64 index = 0; index < edgeCount; ++index) {
            const ui64 rank = ranks[index];
            while (rank >= rowStart + rowLength) {
                rowStart += rowLength;
                ++row;
                rowLength = vertexCount - row;
            }
            const ui64 column = row + (rank - rowStart);

            AppendUInt(&text, verts[row]);
            text.push_back('\t');
            AppendUInt(&text, verts[column]);
            text.push_back('\t');
            AppendUInt(&text, weights[index]);
            text.push_back('\n');
        }

        FILE* output = std::fopen(outputPath.c_str(), "wb");
        if (output == nullptr) {
            throw TGraphCodecError("cannot open output file: " + outputPath);
        }
        if (!text.empty()) {
            std::fwrite(text.data(), 1, text.size(), output);
        }
        std::fclose(output);
    }
}

int main(int argc, char** argv) {
    std::string mode;
    std::string inputPath;
    std::string outputPath;

    try {
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "-s" || argument == "-d") {
                mode = argument;
            } else if (argument == "-i" && index + 1 < argc) {
                inputPath = argv[++index];
            } else if (argument == "-o" && index + 1 < argc) {
                outputPath = argv[++index];
            } else {
                throw TGraphCodecError("unrecognized argument: " + argument);
            }
        }

        if (mode.empty() || inputPath.empty() || outputPath.empty()) {
            throw TGraphCodecError("usage: run -s -i input.tsv -o graph.bin | run -d -i graph.bin -o output.tsv");
        }

        if (mode == "-s") {
            Serialize(inputPath, outputPath);
        } else {
            Deserialize(inputPath, outputPath);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "error: %s\n", error.what());
        return 1;
    }
    return 0;
}
