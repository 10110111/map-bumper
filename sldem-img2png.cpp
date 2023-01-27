#include <iostream>
#include <fstream>
#include <QRegularExpression>
#include <QString>
#include <QImage>
#include <QFile>

constexpr double DEGREE = M_PI/180;
constexpr double LAT_MAX = 60*DEGREE;

double fetch(std::ifstream& in, const off_t width, const off_t height, const off_t rowStride,
             const off_t requestedX, const off_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, off_t(0), height-1);

    using ElemType = float;
    static off_t currentIndexInCacheBegin = 0;
    static off_t currentIndexInCacheEnd = 0;
    static std::vector<ElemType> cache;

    const auto index = x + y*rowStride;
    if(currentIndexInCacheBegin <= index && index < currentIndexInCacheEnd)
        return cache[index - currentIndexInCacheBegin];

    const off_t desiredMaxCacheSize = (1<<30) / sizeof(ElemType);
    const auto numRowsInCache = std::max(off_t(2), desiredMaxCacheSize/rowStride);
    currentIndexInCacheBegin = std::max(off_t(0), index/rowStride - 1) * rowStride;
    currentIndexInCacheEnd = std::min(currentIndexInCacheBegin + rowStride * numRowsInCache,
                                      width*height);
    std::cerr << "Updating cache to range [" << currentIndexInCacheBegin << ", "
              << currentIndexInCacheEnd << ") in numbers, ["
              << currentIndexInCacheBegin/rowStride << ", " << currentIndexInCacheEnd/rowStride
              << ") in rows, total " << numRowsInCache << " rows = "
              << (currentIndexInCacheEnd-currentIndexInCacheBegin)*sizeof(ElemType)
              << " bytes, provoking index: " << index << "\n";

    cache.resize(currentIndexInCacheEnd-currentIndexInCacheBegin);
    in.seekg(currentIndexInCacheBegin * sizeof cache[0]);
    in.read(reinterpret_cast<char*>(cache.data()), cache.size() * sizeof cache[0]);
    return cache[index - currentIndexInCacheBegin];
}

double sample(std::ifstream& in, const off_t width, const off_t height, const off_t rowStride,
              const double longitude, const double latitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-LAT_MAX <= latitude && latitude <= LAT_MAX);

    const auto deltaLon = -2*M_PI/width;
    const auto firstLon = 2*M_PI;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -2*LAT_MAX/(height-1);
    const auto firstLat = LAT_MAX;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetch(in, width, height, rowStride, floorX  , floorY);
    const auto pTopRight    = fetch(in, width, height, rowStride, floorX+1, floorY);
    const auto pBottomLeft  = fetch(in, width, height, rowStride, floorX  , floorY+1);
    const auto pBottomRight = fetch(in, width, height, rowStride, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

int main(int argc, char** argv)
try
{
    if(argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " referenceSphereRadiusInKM kmPerUnitValue inputFile.IMG outputFileName outFileHeight\n";
        return 1;
    }

    const auto outputReferenceRadius = std::stod(argv[1]);
    const auto outputKmPerUnit = std::stod(argv[2]);
    const QString inFileName  = argv[3];
    const char*const outFileName = argv[4];
    const auto outputHeight = std::stoul(argv[5]);

    if(!inFileName.toUpper().endsWith(".IMG"))
    {
        std::cerr << "Input file name must end with .IMG\n";
        return 1;
    }
    const auto labelFileName = QString(inFileName).replace(QRegularExpression("\\.IMG$"), ".LBL");

    QFile file(labelFileName);
    if(!file.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << inFileName.toStdString() << " for reading: "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }
    const auto label = file.readAll();
    if(label.isEmpty())
    {
        std::cerr << "Failed to read label file " << labelFileName.toStdString() << ": "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }
    const auto bytesPerLine = QRegularExpression("RECORD_BYTES\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();
    const auto inputHeight = QRegularExpression("FILE_RECORDS\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();
    const auto inputWidth = bytesPerLine / sizeof(float);
    if(bytesPerLine % sizeof(float))
    {
        std::cerr << "Line size in bytes is not a multiple of sizeof(float)\n";
        return 1;
    }
    if(inputWidth==0 || inputHeight==0)
    {
        std::cerr << "Failed to find input image dimensions\n";
        return 1;
    }

    std::ifstream in(inFileName.toStdString(), std::ios_base::binary);
    if(!in)
    {
        std::cerr << "Failed to open input file\n";
        return 1;
    }
    in.exceptions(std::ios_base::badbit|std::ios_base::failbit);

    const auto outputWidth = 2*outputHeight;
    QImage output(outputWidth, outputHeight, QImage::Format_Grayscale16);
    output.fill(QColor(0,0,0));
    const auto outData = reinterpret_cast<uint16_t*>(output.bits());

    const auto outputRowStride = output.bytesPerLine() / sizeof outData[0];
    const double outputDeltaLat = -M_PI/(outputHeight-1);
    const double outputDeltaLon = -2*M_PI/outputWidth;
    const double outputFirstLat = (1.-outputHeight)/2. * outputDeltaLat;
    const double outputFirstLon = (1.-outputWidth )/2. * outputDeltaLon;
    const double outputLatEpsilon = outputDeltaLat/4; // guaranteed to be larger than rounding error and smaller than step
    for(unsigned j = 0; j < outputHeight; ++j)
    {
        const auto latitude = outputFirstLat + j*outputDeltaLat;
        if(std::abs(latitude) > LAT_MAX + outputLatEpsilon)
            continue;
        for(unsigned i = 0; i < outputWidth; ++i)
        {
            const auto longitude = outputFirstLon + i*outputDeltaLon;
            const auto value = sample(in, inputWidth, inputHeight, inputWidth, longitude, latitude);
            constexpr double inputReferenceRadius = 1737.4; // km
            constexpr double inputKmPerUnit = 1;
            const auto recenteredValue = value + inputReferenceRadius - outputReferenceRadius;
            const auto outputValue = recenteredValue * inputKmPerUnit / outputKmPerUnit;
            outData[j*outputRowStride+i] = outputValue;
        }
    }

    std::cerr << "Saving output file... ";
    if(!output.save(outFileName))
    {
        std::cerr << "Failed to save output file\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
