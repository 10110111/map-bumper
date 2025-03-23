#include <cmath>
#include <thread>
#include <chrono>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <QCryptographicHash>
#include <QImageReader>
#include <QPainter>
#include <QImage>
#include <QDir>
#include "timing.hpp"
#include "healpix.hpp"

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height,
             const size_t rowStride, const size_t channelsPerPixel,
             const size_t subpixelIndex, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);

    const auto rawValue = data[x*channelsPerPixel + subpixelIndex + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height,
              const size_t rowStride, const size_t channelsPerPixel,
              const size_t subpixelIndex, const double longitude, const double latitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

    const auto deltaLon = 2.*M_PI/width;
    const auto firstLon = (1.-width)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -M_PI/height;
    const auto firstLat = (1.-height)/2. * deltaLat;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY+1);

    const auto badLevelLow = 2. / 255.;
    const auto badLevelHigh = 253. / 255.;
    const bool bad = pTopLeft < badLevelLow || pTopRight < badLevelLow || pBottomLeft < badLevelLow || pBottomRight < badLevelLow ||
                     pTopLeft >= badLevelHigh || pTopRight >= badLevelHigh || pBottomLeft >= badLevelHigh || pBottomRight >= badLevelHigh;

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    const auto value = sampleLeft + (sampleRight-sampleLeft)*fracX;
    return bad ? -value : value;
}

template<typename T>
double sampleLROC(T const* data, const size_t width, const size_t height,
                  const size_t rowStride, const size_t yMin, const size_t yMax,
                  const size_t channelsPerPixel, const size_t subpixelIndex,
                  const double longitude, const double latitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

    const auto deltaLon = 2.*M_PI/width;
    const auto firstLon = (1.-width)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -M_PI/(height*90/70);
    const auto firstLat = (1.-height*90/70)/2. * deltaLat;

    auto y = (latitude - firstLat) / deltaLat;
    auto floorY = std::floor(y);

    if(floorY < yMin || floorY+1 >= yMax) return -1;

    y      -= yMin;
    floorY -= yMin;

    const auto pTopLeft     = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

double normalizeLatitude(double x)
{
    if(x < -M_PI/2) return x + 2*M_PI;
    if(x > M_PI/2) return x - 2*M_PI;
    return x;
}

double normalizeLongitude(double x)
{
    if(x < -M_PI) return x + 2*M_PI;
    if(x > M_PI) return x - 2*M_PI;
    return x;
}

struct ColorChannel
{
    enum V
    {
        Red,
        Green,
        Blue,
    } v;
    ColorChannel(V v) : v(v) {}
    ColorChannel& operator=(V v)
    {
        this->v = v;
        return *this;
    }
    operator V() const { return v; }
};

constexpr int HIPS_TILE_SIZE = 512;
const QString initialExt = "tiff"; // Should be uncompressed for saving and reloading speed

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    const auto s = step(0.0031308,c);
    return s  * (1.055*std::pow(c, 1/2.4)-0.055) +
        (1-s) *  12.92*c;
}

double sRGBInverseTransferFunction(const double srgb)
{
	const auto s = step(0.04045, srgb);
	const auto d = 1. - s;
	return s * pow((srgb+0.055)/1.055, float(2.4)) +
	       d * srgb/12.92;
}

double wmsToRed(const double x)
{
    const auto x2 = x*x;
    const auto x4 = x2*x2;
    const auto x8 = x4*x4;
    const auto x16 = x8*x8;
    const auto x32 = x16*x16;
    const auto x48 = x32*x16;
    return 3.55192315861786 + x*(6.56189433868235 + x*(31.7121712083171 + x*(-56.3543748015033 + x*(31.1809216095956 + 1.01475105374168*x48))));
}
double wmsToGreen(const double x)
{
    const auto x2 = x*x;
    const auto x4 = x2*x2;
    const auto x8 = x4*x4;
    const auto x16 = x8*x8;
    const auto x32 = x16*x16;
    const auto x48 = x32*x16;
    return 2.60839445543049 + x*(5.45378170251246 + x*(20.6010742459879 + x*(-35.1688081878114 + x*(19.3469470742881 + 1.00733112382061*x48))));
}
double wmsToBlue(const double x)
{
    const auto x2 = x*x;
    const auto x4 = x2*x2;
    const auto x8 = x4*x4;
    const auto x16 = x8*x8;
    const auto x32 = x16*x16;
    const auto x48 = x32*x16;
    return 1.98418672575485 + x*(4.20391508647678 + x*(13.3557316801473 + x*(-20.5619534557279 + x*(10.9705752833082 + 0.954472577730387*x48))));
}

void fillFace(const int order, const int pix, const uint8_t* data,
              const ssize_t width, const ssize_t height, const ssize_t rowStride,
              const ssize_t channelsPerPixel, const double valueScale,
              const QImage& lrocColorRef, const size_t refYMin, const size_t refYMax,
              const QImage& smallerImg, double* outData)
{
    const unsigned nside = 1u << order;
    int ix, iy, face;
    healpix_nest2xyf(nside, pix, &ix, &iy, &face);

    const auto& ref = lrocColorRef;
    const auto& small = smallerImg;

    for(int y = 0; y < HIPS_TILE_SIZE; ++y)
    {
        for(int x = 0; x < HIPS_TILE_SIZE; ++x)
        {
            double theta, phi;
            healpix_xyf2ang(nside * HIPS_TILE_SIZE,
                            ix * HIPS_TILE_SIZE + x, iy * HIPS_TILE_SIZE + y,
                            face, &theta, &phi);
            const double latitude = normalizeLatitude(M_PI/2 - theta);
            const double longitude = normalizeLongitude(M_PI+phi);
            assert(-M_PI <= longitude && longitude <= M_PI);
            assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

            // HiPS coordinates are swapped, because they were designed from the "look from
            // sphere center" perspective, while we are making a look from outside a planet.
            const int i = y, j = x;

            const int pixelPosInOutData = (i + j*HIPS_TILE_SIZE)*channelsPerPixel;
            const auto samp = sample(data, width, height, rowStride, 1, 0, longitude, latitude);
            const auto sampSm = sample(small.bits(), small.width(), small.height(), small.bytesPerLine(), 1, 0, longitude, latitude);
            const auto refRsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 0, longitude, latitude);
            const auto refGsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 1, longitude, latitude);
            const auto refBsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 2, longitude, latitude);
            const auto refR = sRGBInverseTransferFunction(refRsRGB);
            const auto refG = sRGBInverseTransferFunction(refGsRGB);
            const auto refB = sRGBInverseTransferFunction(refBsRGB);
            const double smallR = wmsToRed(std::abs(sampSm));
            const double smallG = wmsToGreen(std::abs(sampSm));
            const double smallB = wmsToBlue(std::abs(sampSm));
            const double red   = wmsToRed(std::abs(samp));
            const double green = wmsToGreen(std::abs(samp));
            const double blue  = wmsToBlue(std::abs(samp));
            if(refRsRGB < 0 || refGsRGB < 0 || refBsRGB < 0)
            {
                outData[pixelPosInOutData + 0] = sRGBTransferFunction(std::clamp(red   * valueScale, 0., 1.));
                outData[pixelPosInOutData + 1] = sRGBTransferFunction(std::clamp(green * valueScale, 0., 1.));
                outData[pixelPosInOutData + 2] = sRGBTransferFunction(std::clamp(blue  * valueScale, 0., 1.));
            }
            else if(samp < 0)
            {
                outData[pixelPosInOutData + 0] = refRsRGB;
                outData[pixelPosInOutData + 1] = refGsRGB;
                outData[pixelPosInOutData + 2] = refBsRGB;
            }
            else
            {
                const auto R = red, G = green, B = blue;
                const auto refI = refR+refG+refB;
                const auto smallI = smallR+smallG+smallB;
                const auto I = (R+G+B) * refI / smallI;
                const auto rRef = refR/refI, gRef = refG/refI, bRef = refB/refI;
                // TODO: adjust I so that the average of the current 4×4 cell is the
                // same as that of the corresponding pixel of the reference image.
                outData[pixelPosInOutData + 0] = sRGBTransferFunction(std::clamp(rRef * I, 0., 1.));
                outData[pixelPosInOutData + 1] = sRGBTransferFunction(std::clamp(gRef * I, 0., 1.));
                outData[pixelPosInOutData + 2] = sRGBTransferFunction(std::clamp(bRef * I, 0., 1.));
            }
        }
    }
}

void createLowerOrderTile(const int order, const int pix, const QString& outDir, QImage::Format imgFormat)
{
    try
    {
        const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(order).arg((pix / 10000) * 10000);
        if(!QDir().mkpath(outPath))
            throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
        const auto outFileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(initialExt);
        const auto inFileTemplate = QString("%1/Norder%2/Dir%4/Npix%5.%3").arg(outDir).arg(order + 1).arg(initialExt);

        QImage outImg(HIPS_TILE_SIZE, HIPS_TILE_SIZE, imgFormat);
        {
            QPainter p(&outImg);
            for(int j = 0; j < 2; ++j)
            {
                for(int i = 0; i < 2; ++i)
                {
                    const int innerPix = pix * 4 + i * 2 + j;
                    const auto path = inFileTemplate.arg((innerPix / 10000) * 10000).arg(innerPix);
                    QImage img(path);
                    if(img.isNull())
                        throw std::runtime_error("Failed to open \""+path.toStdString()+'"');
                    img = img.scaled(HIPS_TILE_SIZE / 2, HIPS_TILE_SIZE / 2, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
                    p.drawImage(QPoint(HIPS_TILE_SIZE / 2 * i, HIPS_TILE_SIZE / 2 * j), img);
                }
            }
        }
        outImg.save(outFileName, nullptr, 100);
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << "\n";
    }
}

QString md5sum(QString const& filePath)
{
    QFile file(filePath);
    if(!file.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open \"" << filePath.toStdString() << "\" for reading to compute MD5 sum\n";
        goto error;
    }
    {
        QCryptographicHash hash(QCryptographicHash::Md5);
        if(hash.addData(&file))
            return hash.result().toHex();
    }
    std::cerr << "Failed to compute MD5 sum of \"" << filePath.toStdString() << "\"\n";
error:
    std::cerr << "The property \"source_md5\" will be left empty\n";
    return "";
}

QString formatDate()
{
    return QDateTime::currentDateTime().toUTC().toString("yyyy-MM-ddTHH:mmZ");
}

void handleProgressReporting(const size_t totalItemCount,
                             const std::chrono::time_point<std::chrono::steady_clock> startTime,
                             std::chrono::time_point<std::chrono::steady_clock>& time0,
                             std::atomic_int& numThreadsReportedFirstProgress,
                             size_t& itemsDoneInThisThreadAfterLastUpdate,
                             std::atomic<unsigned>& itemsDone)
{
    ++itemsDoneInThisThreadAfterLastUpdate;
    using namespace std::chrono;
    const auto time1 = steady_clock::now();
    if(time1 - time0 > seconds(5))
    {
        // ETA is reported only by the first thread in the iteration.
        // This seems to provide the most accurate estimate. Once the
        // reporting thread is determined, it will report after each
        // iteration.
        thread_local bool isTimeReportingThread = false;
        thread_local bool reportedProgressFirstTime = false;
        if(!reportedProgressFirstTime)
        {
            if(numThreadsReportedFirstProgress.fetch_add(1) == 0)
                isTimeReportingThread = true;
        }

        itemsDone += itemsDoneInThisThreadAfterLastUpdate;
        itemsDoneInThisThreadAfterLastUpdate = 0;
        const auto progress = std::round(itemsDone * 10000. / totalItemCount)/10000.;
        const auto usecElapsed = duration_cast<microseconds>(time1 - startTime).count();
        const long secToEnd = std::lround((1. - progress) / progress * usecElapsed * 1e-6);
        std::ostringstream ss;
        ss << progress*100 << "% done";
        // First ETA estimate is too far from reality, likely due
        // to overhead of thread startup, so skip first measurement.
        if(isTimeReportingThread && reportedProgressFirstTime)
        {
            if(secToEnd > 60)
            {
                const auto min = secToEnd/60;
                const auto sec = secToEnd - min*60;
                ss << ", ETA: " << min << 'm' << sec << "s\n";
            }
            else
            {
                ss << ", ETA: " << secToEnd << "s\n";
            }
        }
        else
        {
            ss << '\n';
        }
        reportedProgressFirstTime = true;
        std::cerr << ss.str();
        time0 = time1;
    }
}

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} inDir outDir";
    s << R"(
Options:
 -h, --help                 This help message
 --value-scale N            Scale albedo by N before saving as color
 --lroc-image filename.ext  Use a LROC-derived sRGB 27360×10640 image filename.ext to correct colors
 --small-data filename.ext  Use a reduced-resolution version of the input map for correction
                             of brightness (used together with the --lroc-image map)
 --format FORMAT            Set hips_tile_format to FORMAT (only one value is supported)
 --title TITLE              Set obs_title to TITLE
 --type TYPE                Set type (a non-standard property) to TYPE
 -d, --desc DESCRIPTION     Set obs_description to DESCRIPTION
 --frame FRAME              Set hips_frame to FRAME
 --creator CREATOR          Set hips_creator to CREATOR
 --my-copyright COPYRIGHT   Set hips_copyright to COPYRIGHT
 --orig-copyright COPYRIGHT Set obs_copyright to COPYRIGHT
 -s, --status STATUS        Set hips_status to STATUS
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    setenv("QT_IMAGEIO_MAXALLOC","2048", false);

    QString inDir;
    QString outDir;
    QString imgFormat = "png";
    QString surveyTitle;
    QString surveyType;
    QString description;
    QString frame;
    QString creator;
    QString hips_copyright;
    QString obs_copyright;
    QString hipsStatus = "public mirror clonable";
    ColorChannel colorChannel = ColorChannel::Green;
    constexpr double defaultScale = 105;
    double valueScale = 1 / defaultScale;
    QString lrocSRGBFileName;
    QString smallerImgFileName;

    int totalPositionalArgumentsFound = 0;
    for(int n = 1; n < argc; ++n)
    {
        if(argv[n][0]!='-')
        {
            // Must be a positional argument
            switch(totalPositionalArgumentsFound)
            {
            case 0:
                inDir = argv[n];
                break;
            case 1:
                outDir = argv[n];
                break;
            default:
                std::cerr << "Extraneous positional argument\n";
                return usage(argv[0], 1);
            }
            ++totalPositionalArgumentsFound;
            continue;
        }
        // OK, we got a switch
        const std::string arg = argv[n];
        if(arg == "-h" || arg == "--help")
            return usage(argv[0], 0);

#define GO_TO_PARAM()                                                       \
            ++n;                                                            \
            if(n == argc)                                                   \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return 1;                                                   \
            }                                                               \
            do{}while(0)

        if(arg == "--format")
        {
            GO_TO_PARAM();
            imgFormat = argv[n];
        }
        else if(arg == "--title")
        {
            GO_TO_PARAM();
            surveyTitle = argv[n];
        }
        else if(arg == "--type")
        {
            GO_TO_PARAM();
            surveyType = argv[n];
        }
        else if(arg == "-d" || arg == "--desc" || arg == "--description")
        {
            GO_TO_PARAM();
            description = argv[n];
        }
        else if(arg == "--frame")
        {
            GO_TO_PARAM();
            frame = argv[n];
        }
        else if(arg == "--creator")
        {
            GO_TO_PARAM();
            creator = argv[n];
        }
        else if(arg == "--my-copyright")
        {
            GO_TO_PARAM();
            hips_copyright = argv[n];
        }
        else if(arg == "--orig-copyright")
        {
            GO_TO_PARAM();
            obs_copyright = argv[n];
        }
        else if(arg == "-s" || arg == "--status")
        {
            GO_TO_PARAM();
            hipsStatus = argv[n];
        }
        else if(arg == "--value-scale")
        {
            GO_TO_PARAM();
            valueScale = std::stod(argv[n]) / defaultScale;
        }
        else if(arg == "--lroc-image")
        {
            GO_TO_PARAM();
            lrocSRGBFileName = argv[n];
        }
        else if(arg == "--small-data")
        {
            GO_TO_PARAM();
            smallerImgFileName = argv[n];
        }
        else if(arg == "--red")
        {
            colorChannel = ColorChannel::Red;
        }
        else if(arg == "--green")
        {
            colorChannel = ColorChannel::Green;
        }
        else if(arg == "--blue")
        {
            colorChannel = ColorChannel::Blue;
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(inDir.isEmpty())
    {
        std::cerr << "Input file not specified\n";
        return usage(argv[0], 1);
    }
    if(outDir.isEmpty())
    {
        std::cerr << "Output directory not specified\n";
        return usage(argv[0], 1);
    }
    if(lrocSRGBFileName.isEmpty())
    {
        std::cerr << "LROC-derived sRGB image not specified\n";
        return 1;
    }
    if(smallerImgFileName.isEmpty())
    {
        std::cerr << "Reduced-resolution map not specified\n";
        return 1;
    }

    QString finalExt;
    if(imgFormat == "jpeg")
        finalExt = "jpg";
    else if(imgFormat == "webp" || imgFormat == "bmp" || imgFormat == "png" || imgFormat == "tiff")
        finalExt = imgFormat;
    else
        throw std::runtime_error("Unexpected output image format: " + imgFormat.toStdString());

    std::cerr << "Preparing data container...";
    constexpr ssize_t inputTileWidth = 4962;
    constexpr ssize_t inputTileHeight = inputTileWidth;
    constexpr ssize_t tileCountVert = 11;
    constexpr ssize_t tileCountHoriz = tileCountVert * 2;
    constexpr ssize_t height = inputTileHeight * tileCountVert;
    constexpr ssize_t width = inputTileWidth * tileCountHoriz;
    std::unique_ptr<uint8_t[]> bigData(new uint8_t[height*width]);
    std::cerr << "done\n";

    std::cerr << "Loading the reduced-resolution map... ";
    QImage smallerImg(smallerImgFileName);
    if(smallerImg.isNull())
        throw std::runtime_error(("Failed to load image \""+smallerImgFileName+"\"").toStdString());
    smallerImg = smallerImg.convertToFormat(QImage::Format_Grayscale8);
    if(smallerImg.isNull())
        throw std::runtime_error("Failed to convert the reduced-resolution map to the working format");
    std::cerr << "done\n";

    std::cerr << "Loading LROC-derived sRGB image... ";
    QImage lrocImg(lrocSRGBFileName);
    if(lrocImg.isNull())
        throw std::runtime_error(("Failed to load image \""+lrocSRGBFileName+"\"").toStdString());
    lrocImg = lrocImg.convertToFormat(QImage::Format_RGB888);
    if(lrocImg.isNull())
        throw std::runtime_error("Failed to convert the LROC-derived sRGB image to the working format");
    std::cerr << "done\n";

    std::cerr << "Loading input images...\n";
    for(int tileJ = 0; tileJ < tileCountVert; ++tileJ)
    {
        for(int tileI = 0; tileI < tileCountHoriz; ++tileI)
        {
            const auto inFileName = QString("tile-%1-%2.tiff").arg(tileI,2,10,QChar('0')).arg(tileJ,2,10,QChar('0'));
            std::cerr << " " << inFileName.toStdString() << " ";
            const auto inFilePath = inDir + "/" + inFileName;
            QImage img(inFilePath);
            if(img.isNull())
                throw std::runtime_error(("Failed to open input image \""+inFilePath+"\"").toStdString());
            if(img.width() != inputTileWidth || img.height() != inputTileWidth)
            {
                throw std::runtime_error("Unexpected tile dimensions: expected "+
                                         std::to_string(inputTileWidth)+u8"×"+std::to_string(inputTileWidth)+
                                         ", got "+std::to_string(img.width())+u8"×"+std::to_string(img.height()));
            }
            img = img.convertToFormat(QImage::Format_Grayscale8);
            if(img.isNull())
                throw std::runtime_error("Failed to convert input image to the working format");
            const auto tilePos = inputTileHeight*(tileCountVert - tileJ - 1)*width + inputTileWidth*tileI;
            for(int bigImgJ = 0; bigImgJ < inputTileHeight; ++bigImgJ)
            {
                std::memcpy(bigData.get() + tilePos + bigImgJ*width,
                            img.bits() + bigImgJ*img.bytesPerLine(),
                            inputTileWidth);
            }
            std::cerr << "\n";
        }
    }
    std::cerr << "done\n";

    const int orderMax = std::ceil(std::log2(width / (4. * HIPS_TILE_SIZE * M_SQRT2)));
    // First create the tiles of the deepest level
    std::cerr << "Creating tiles of order " << orderMax << "...\n";
    {
        const int absolutePixMax = 12 * (1 << (2 * orderMax));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,orderMax,outDir,startTime,lrocImg,smallerImg,
                     &bigData = std::as_const(bigData),colorChannel, valueScale,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            constexpr int channelsPerPixel = 3;
            std::vector<double> data(HIPS_TILE_SIZE * HIPS_TILE_SIZE * channelsPerPixel);
            const size_t refYmin = (lrocImg.height()*90/70 - lrocImg.height()) / 2;
            const size_t refYmax = lrocImg.height()*90/70 - refYmin;

            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                fillFace(orderMax, pix, bigData.get(), width, height, width, channelsPerPixel,
                         valueScale, lrocImg, refYmin, refYmax, smallerImg, data.data());

                using OutType = uchar;
                std::vector<OutType> outBits;
                outBits.reserve(data.size());
                if(std::is_integral_v<OutType>)
                {
                    for(auto v : data)
                        outBits.push_back(std::clamp(v,0.,1.)*std::numeric_limits<OutType>::max());
                }
                else
                {
                    for(auto v : data)
                        outBits.push_back(v*std::numeric_limits<OutType>::max());
                }
                const QImage out(outBits.data(), HIPS_TILE_SIZE, HIPS_TILE_SIZE, channelsPerPixel*HIPS_TILE_SIZE, QImage::Format_RGB888);
                const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(orderMax).arg((pix / 10000) * 10000);
                if(!QDir().mkpath(outPath))
                    throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
                const auto fileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(initialExt);
                if(!out.save(fileName, nullptr, 100))
                    throw std::runtime_error("Failed to save output file " + fileName.toStdString());

                handleProgressReporting(absolutePixMax, startTime, time0, numThreadsReportedFirstProgress,
                                        itemsDoneInThisThreadAfterLastUpdate, itemsDone);
            }
        };
        const auto time0 = std::chrono::steady_clock::now();
        const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        for(size_t n = 0; n < numThreads; ++n)
        {
            const size_t pixMin = absolutePixMax / numThreads * n;
            const size_t pixMax = n+1 < numThreads ? absolutePixMax / numThreads * (n+1) : absolutePixMax;
            threads.emplace_back(work, pixMin, pixMax);
        }
        for(auto& thread : threads)
            thread.join();
        auto time1 = std::chrono::steady_clock::now();
        std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";
    }

    // Now merge the deepest-level tiles to create the ones with lower detail level
    for(int order = orderMax - 1; order >= 0; --order)
    {
        std::cerr << "Creating tiles of order " << order << "...\n";
        const int absolutePixMax = 12 * (1 << (2 * order));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,order,outDir,startTime,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                createLowerOrderTile(order, pix, outDir, QImage::Format_RGB888);
                handleProgressReporting(absolutePixMax, startTime, time0, numThreadsReportedFirstProgress,
                                        itemsDoneInThisThreadAfterLastUpdate, itemsDone);
            }
        };
        const auto time0 = std::chrono::steady_clock::now();
        const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        for(size_t n = 0; n < numThreads; ++n)
        {
            const size_t pixMin = absolutePixMax / numThreads * n;
            const size_t pixMax = n+1 < numThreads ? absolutePixMax / numThreads * (n+1) : absolutePixMax;
            threads.emplace_back(work, pixMin, pixMax);
        }
        for(auto& thread : threads)
            thread.join();
        auto time1 = std::chrono::steady_clock::now();
        std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";
    }

    if(finalExt != initialExt)
    {
        std::cerr << "Converting tiles to " << imgFormat.toStdString() << "...\n";
        struct OrderAndPix
        {
            int order;
            int pix;
            OrderAndPix(int order, int pix) : order(order), pix(pix) {}
        };
        std::vector<OrderAndPix> jobsToDo;
        for(int order = 0; order <= orderMax; ++order)
        {
            const int pixMax = 12 * (1 << (2 * order));
            for(int pix = 0; pix < pixMax; ++pix)
                jobsToDo.emplace_back(order, pix);
        }
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [&jobsToDo=std::as_const(jobsToDo),outDir,startTime,finalExt,
                     &numThreadsReportedFirstProgress,&itemsDone](const int jobMin, const int jobMax)
         {
             auto time0 = std::chrono::steady_clock::now();
             size_t itemsDoneInThisThreadAfterLastUpdate = 0;
             for(int n = jobMin; n < jobMax; ++n)
             {
                 const int order = jobsToDo[n].order;
                 const int pix = jobsToDo[n].pix;
                 const auto pathTemplate = QString("%1/Norder%2/Dir%3/Npix%4.%5").arg(outDir).arg(order).arg((pix / 10000) * 10000).arg(pix);
                 const auto inFileName = pathTemplate.arg(initialExt);
                 const auto outFileName = pathTemplate.arg(finalExt);
                 QImage img(inFileName);
                 if(img.isNull())
                 {
                     std::cerr << "Failed to open recently saved file " << inFileName.toStdString() << "\n";
                     break;
                 }
                 if(!img.save(outFileName))
                 {
                     std::cerr << "Failed to save output file "  << outFileName.toStdString() << "\n";
                     break;
                 }
                 if(!QFile(inFileName).remove())
                     std::cerr << "Warning: failed to remove " << inFileName.toStdString() << "\n";

                 handleProgressReporting(jobsToDo.size(), startTime, time0, numThreadsReportedFirstProgress,
                                         itemsDoneInThisThreadAfterLastUpdate, itemsDone);
             }
         };
        const auto time0 = std::chrono::steady_clock::now();
        const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        for(size_t n = 0; n < numThreads; ++n)
        {
            const size_t jobMin = jobsToDo.size() / numThreads * n;
            const size_t jobMax = n+1 < numThreads ? jobsToDo.size() / numThreads * (n+1) : jobsToDo.size();
            threads.emplace_back(work, jobMin, jobMax);
        }
        for(auto& thread : threads)
            thread.join();
        auto time1 = std::chrono::steady_clock::now();
        std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";
    }

    const auto propsPath = outDir+"/properties";
    QFile propsFile(propsPath);
    if(!propsFile.open(QFile::WriteOnly))
        throw std::runtime_error("Failed to open \""+propsPath.toStdString()+"\" for writing");
    {
        QTextStream props(&propsFile);
        props << "hips_order            = " << orderMax << "\n";
        props << "hips_order_min        = 0\n";
        props << "hips_tile_width       = " << HIPS_TILE_SIZE << "\n";
        props << "hips_tile_format      = " << imgFormat << "\n";
        props << "dataproduct_type      = image\n";
        props << "obs_title             = " << surveyTitle << "\n";
        props << "hips_release_date     = " << formatDate() << "\n";
        if(!surveyType.isEmpty())
            props << "type                  = " << surveyType << "\n";
        if(!description.isEmpty())
            props << "obs_description       = " << description << "\n";
        props << "hips_frame            = " << frame << "\n";
        if(!obs_copyright.isEmpty())
            props << "obs_copyright         = " << obs_copyright << "\n";
        if(!hips_copyright.isEmpty())
            props << "hips_copyright        = " << hips_copyright << "\n";
        if(!creator.isEmpty())
            props << "hips_creator          = " << creator << "\n";
        props << "hips_version          = 1.4\n";
        props << "hips_status           = " << hipsStatus << "\n";
    }
    if(!propsFile.flush())
        throw std::runtime_error("Failed to write properties file");

    if(frame.isEmpty())
        std::cerr << "Warning: hips_frame is not specified\n";
    if(surveyTitle.isEmpty())
        std::cerr << "Warning: obs_title is not specified\n";
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
