#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <QRegularExpression>
#include <QString>
#include <QImage>
#include <QFile>

constexpr double DEGREE = M_PI/180;

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
              const double longitude, const double latitude, const double maxLatitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-maxLatitude <= latitude && latitude <= maxLatitude);

    const auto deltaLon = -2*M_PI/width;
    const auto firstLon = 2*M_PI;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -2*maxLatitude/(height-1);
    const auto firstLat = maxLatitude;

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

int usage(const char*const argv0, const int ret)
{
        std::cerr << "Usage: " << argv0 << 1+R"(
 {-h outputMapHeight|-w outputMapWidth|-b backgroundAltitudeMap} -r referenceSphereRadiusInKM -k kmPerUnitValue -i inputFile.IMG -o outputFileName
Options:
 --help                 Show this help message and quit
 -b,--background FILE   Background altitude map, a map that has full -90°..90° latitude span, in which the -60°..60° region is to be replaced by the data from input file
 -h,--height NUM        Height of the output map (width is 2×this)
 -w,--width NUM         Width of the output map (must be even; height is 0.5×this)
 -r,--radius NUM        Radius of the reference sphere in km
 -k,--km-per-unit NUM   Unit of altitude to use in the output
 -l,--max-lat NUM       Override maximum latitude (in degrees) by rescaling the map (default is 60° as for SLDEM2015)
 -i FILE                Input file from SLDEM2015, must be accompanied by the corresponding .LBL file
 -o FILE                Output file to save the resulting map in
)";
    return ret;
}


int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    double outputReferenceRadius = -1;
    double outputKmPerUnit = -1;
    off_t outputWidth = 0, outputHeight = 0;
    QString inFileName;
    QString outFileName;
    QString bgMapFileName;
    double maxLatitude = 60*DEGREE;

    if(argc == 1)
        return usage(argv[0],1);

    for(int n = 1; n < argc; ++n)
    {
#define REQUIRE_PARAM() do{                                                 \
            if(argc+1 == n)                                                 \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return usage(argv[0], 1);                                   \
            }}while(false)

        const std::string arg = argv[n];
        if(arg == "--help")
        {
            return usage(argv[0], 0);
        }
        else if(arg == "-b" || arg == "--background")
        {
            REQUIRE_PARAM();
            bgMapFileName = argv[++n];
        }
        else if(arg == "-h" || arg == "--height")
        {
            REQUIRE_PARAM();
            outputHeight = std::stoul(argv[++n]);
        }
        else if(arg == "-w" || arg == "--width")
        {
            REQUIRE_PARAM();
            outputWidth = std::stoul(argv[++n]);
        }
        else if(arg == "-r" || arg == "--radius")
        {
            REQUIRE_PARAM();
            outputReferenceRadius = std::stod(argv[++n]);
        }
        else if(arg == "-k" || arg == "--km-per-unit")
        {
            REQUIRE_PARAM();
            outputKmPerUnit = std::stod(argv[++n]);
        }
        else if(arg == "-l" || arg == "--max-lat")
        {
            REQUIRE_PARAM();
            maxLatitude = std::stod(argv[++n])*DEGREE;
        }
        else if(arg == "-i")
        {
            REQUIRE_PARAM();
            inFileName = argv[++n];
        }
        else if(arg == "-o")
        {
            REQUIRE_PARAM();
            outFileName = argv[++n];
        }
        else
        {
            std::cerr << "Unknown option " << arg << "\n";
            return usage(argv[0], 1);
        }
    }
#undef REQUIRE_PARAM

    if(outputWidth <= 0 && outputHeight <= 0 && bgMapFileName.isEmpty())
    {
        std::cerr << "Neither output width nor height, nor background altitude map were specified\n";
        return usage(argv[0],1);
    }
    if(outputHeight > 0 && outputWidth > 0 && outputWidth != outputHeight*2)
    {
        std::cerr << "Width must be twice the height. If one is specified the other may be omitted.\n";
        return 1;
    }
    if(outputWidth % 2)
    {
        std::cerr << "Width must be even\n";
        return 1;
    }
    if(outputReferenceRadius <= 0)
    {
        std::cerr << "Reference radius was not specified\n";
        return usage(argv[0],1);
    }
    if(outputKmPerUnit <= 0)
    {
        std::cerr << "Altitude unit was not specified\n";
        return usage(argv[0],1);
    }
    if(inFileName.isEmpty() || outFileName.isEmpty())
    {
        std::cerr << (inFileName.isEmpty() ? "Input" : "Output") << " file was not specified\n";
        return usage(argv[0],1);
    }

    if(!inFileName.toUpper().endsWith(".IMG"))
    {
        std::cerr << "Input file name must end with .IMG\n";
        return 1;
    }
    const auto labelFileName = inFileName.endsWith(".IMG") ? QString(inFileName).replace(QRegularExpression("\\.IMG$"), ".LBL")
                                                           : QString(inFileName).replace(QRegularExpression("\\.img$"), ".lbl");

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

    QImage output;
    if(bgMapFileName.isEmpty())
    {
        if(outputWidth <= 0)
            outputWidth = outputHeight * 2;
        if(outputHeight <= 0)
            outputHeight = outputWidth / 2;

        output = QImage(outputWidth, outputHeight, QImage::Format_Grayscale16);
        output.fill(QColor(0,0,0));
    }
    else
    {
        output = QImage(bgMapFileName).convertToFormat(QImage::Format_Grayscale16);
        if(output.isNull())
        {
            std::cerr << "Failed to load background altitude map\n";
            return 1;
        }
        if(outputWidth && outputWidth != output.width())
        {
            std::cerr << "Output width doesn't match background altitude map width\n";
            return 1;
        }
        if(outputHeight && outputHeight != output.height())
        {
            std::cerr << "Output height doesn't match background altitude map height\n";
            return 1;
        }
        // For the case they were not specified
        outputWidth = output.width();
        outputHeight = output.height();
    }
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
        if(std::abs(latitude) > maxLatitude + outputLatEpsilon)
            continue;
        for(unsigned i = 0; i < outputWidth; ++i)
        {
            const auto longitude = outputFirstLon + i*outputDeltaLon;
            const auto value = sample(in, inputWidth, inputHeight, inputWidth, longitude, latitude, maxLatitude);
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
