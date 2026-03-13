#include <cmath>
#include <thread>
#include <chrono>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <QImageReader>
#include <QImage>
#include <QDir>
#include "hips.hpp"
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

double empToRed(const double x, const double lon, const double lat)
{
    constexpr double PI = M_PI;
    if(lat < 0)
    {
        const double r = (PI/2 + lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (9.08617789160424 + r*cos(lon) + 0.2094270102229348*r2*cos(2*(lon + PI)) + 0.1490437731690924*r3*cos(3*(lon + PI)) -
                0.2789006239867136*r*sin(lon) - 0.5014459534312736*r2*sin(2*(lon + PI)) - 0.5129670453772869*r3*sin(3*(lon + PI)))
                                                                    /
               (5.919064090108278 +
                0.65194174632201190*r*cos(lon) + 0.20942016254585352*r2*cos(2*(lon + PI)) + 0.10057140731911288*r3*cos(3*(lon + PI)) -
                0.16336596272597545*r*sin(lon) - 0.27991274292544255*r2*sin(2*(lon + PI)) - 0.33717656278579110*r3*sin(3*(lon + PI)))
                                                                    ;
    }
    else
    {
        const double r = (PI/2 - lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (18.034454013991404 + r*cos(lon) - 0.4345698349619620*r2*cos(2*(lon + PI)) - 0.1730334474676742*r3*cos(3*(lon + PI)) +
                0.6830217975120059 * r*sin(lon) - 0.6586493467037962*r2*sin(2*(lon + PI)) + 0.8273972743674709*r3*sin(3*(lon + PI)))
                                                                    /
               (11.864375187315154 +
                0.76541540001983150*r*cos(lon) - 0.23153832556531215*r2*cos(2*(lon + PI)) - 0.12599561590820352*r3*cos(3*(lon + PI)) +
                0.48506666131023424*r*sin(lon) - 0.46478065132830876*r2*sin(2*(lon + PI)) + 0.50390690546784380*r3*sin(3*(lon + PI)))
                                                                    ;
    }
}
double empToGreen(const double x, const double lon, const double lat)
{
    constexpr double PI = M_PI;
    if(lat < 0)
    {
        const double r = (PI/2 + lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (8.314560480826950+r*cos(lon) + 0.19178740432702773*r2*cos(2*(lon + PI)) + 0.16971768410906835*r3*cos(3*(lon + PI)) -
                0.292226486858174*r*sin(lon) - 0.51746044230973050*r2*sin(2*(lon + PI)) - 0.48463412059223643*r3*sin(3*(lon + PI)))
                                                                    /
               (7.203142849483428 +
                0.79337365786374410*r*cos(lon) + 0.25485166631338624*r2*cos(2*(lon + PI)) + 0.12238931737600109*r3*cos(3*(lon + PI)) -
                0.19880649176026408*r*sin(lon) - 0.34063687130068077*r2*sin(2*(lon + PI)) - 0.41032347517619840*r3*sin(3*(lon + PI)))
                                                                    ;
    }
    else
    {
        const double r = (PI/2 - lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (16.299511190548507+r*cos(lon) - 0.40891051668851935*r2*cos(2*(lon + PI)) - 0.18106562828383546*r3*cos(3*(lon + PI)) +
                0.7052590279928693*r*sin(lon) - 0.62250920969539040*r2*sin(2*(lon + PI)) + 0.78809756405954090*r3*sin(3*(lon + PI)))
                                                                    /
               (14.2274852454926 +
                0.9178685045376844*r*cos(lon) - 0.2776554229563270*r2*cos(2*(lon + PI)) - 0.15109103834201776*r3*cos(3*(lon + PI)) +
                0.5816807592405072*r*sin(lon) - 0.5573542436717532*r2*sin(2*(lon + PI)) + 0.60427354575828810*r3*sin(3*(lon + PI)))
                                                                    ;
    }
}
double empToBlue(const double x, const double lon, const double lat)
{
    constexpr double PI = M_PI;
    if(lat < 0)
    {
        const double r = (PI/2 + lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (7.5182390052249750+r*cos(lon) + 0.16416593369496188*r2*cos(2*(lon + PI)) + 0.18193794123524465*r3*cos(3*(lon + PI)) -
                0.2740888215065235*r*sin(lon) - 0.51829680361104250*r2*sin(2*(lon + PI)) - 0.45091205800768436*r3*sin(3*(lon + PI)))
                                                                    /
               (8.683697038188974 +
                0.95644590520131360*r*cos(lon) + 0.30723459275860193*r2*cos(2*(lon + PI)) + 0.14754556101579588*r3*cos(3*(lon + PI)) -
                0.23966973580083203*r*sin(lon) - 0.41065232943753310*r2*sin(2*(lon + PI)) - 0.49466251336977020*r3*sin(3*(lon + PI)))
                                                                    ;
    }
    else
    {
        const double r = (PI/2 - lat) / ((90 - 70) * PI / 180);
        const double r2 = r*r;
        const double r3 = r2*r;
        return x *
               (14.026755652524963+r*cos(lon) - 0.3884757975572265*r2*cos(2*(lon + PI)) - 0.16588908530583166*r3*cos(3*(lon + PI)) +
                0.6997034457568012*r*sin(lon) - 0.5853193082785430*r2*sin(2*(lon + PI)) + 0.73056993587859400*r3*sin(3*(lon + PI)))
                                                                    /
               (16.329657978362214 +
                1.0534875622492514*r*cos(lon) - 0.31868021751424686*r2*cos(2*(lon + PI)) - 0.17341539542291312*r3*cos(3*(lon + PI)) +
                0.6676266175711410*r*sin(lon) - 0.63970575368067180*r2*sin(2*(lon + PI)) + 0.69355758641402880*r3*sin(3*(lon + PI)))
                                                                    ;
    }
}

constexpr int DOWNSAMP_FACTOR = 4;

void fillPoint(const double longitude, const double latitude, const int pixelPosInOutData,
               const uint8_t*const data, const ssize_t width, const ssize_t height, const ssize_t rowStride,
               const QImage& lrocColorRef, const size_t refYMin, const size_t refYMax,
               const uint8_t*const smallData, double* outData)
{
    const auto& ref = lrocColorRef;

    const auto samp = sRGBInverseTransferFunction(sample(data, width, height, rowStride, 1, 0, longitude, latitude));
    const auto sampSm = sRGBInverseTransferFunction(sample(smallData, width/DOWNSAMP_FACTOR, height/DOWNSAMP_FACTOR, rowStride/DOWNSAMP_FACTOR,
                                                           1, 0, longitude, latitude));
    const auto refRsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 0, longitude, latitude);
    const auto refGsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 1, longitude, latitude);
    const auto refBsRGB = sampleLROC(ref.bits(), ref.width(), ref.height(), ref.bytesPerLine(), refYMin, refYMax, 3, 2, longitude, latitude);
    const auto refR = sRGBInverseTransferFunction(refRsRGB);
    const auto refG = sRGBInverseTransferFunction(refGsRGB);
    const auto refB = sRGBInverseTransferFunction(refBsRGB);
    const double smallR = empToRed(std::abs(sampSm), longitude, latitude);
    const double smallG = empToGreen(std::abs(sampSm), longitude, latitude);
    const double smallB = empToBlue(std::abs(sampSm), longitude, latitude);
    const double red   = empToRed(std::abs(samp), longitude, latitude);
    const double green = empToGreen(std::abs(samp), longitude, latitude);
    const double blue  = empToBlue(std::abs(samp), longitude, latitude);
    const double fallbackRed = 0.702, fallbackGreen = 0.616, fallbackBlue = 0.541;
    if(refRsRGB < 0 || refGsRGB < 0 || refBsRGB < 0)
    {
        // No Hapke-normalized data, use a colorized empirically-normalized point
        if(samp > 0)
        {
            outData[pixelPosInOutData + 0] = sRGBTransferFunction(std::clamp(red  , 0., 1.));
            outData[pixelPosInOutData + 1] = sRGBTransferFunction(std::clamp(green, 0., 1.));
            outData[pixelPosInOutData + 2] = sRGBTransferFunction(std::clamp(blue , 0., 1.));
        }
        else
        {
            // No empirically-normalized data either (happens in many craters near the poles)
            outData[pixelPosInOutData + 0] = fallbackRed;
            outData[pixelPosInOutData + 1] = fallbackGreen;
            outData[pixelPosInOutData + 2] = fallbackBlue;
        }
    }
    else if(samp < 0)
    {
        // No empirically-normalized data, use the Hapke-normalized color
        outData[pixelPosInOutData + 0] = refRsRGB;
        outData[pixelPosInOutData + 1] = refGsRGB;
        outData[pixelPosInOutData + 2] = refBsRGB;
    }
    else
    {
        // Normal case: both empirically-normalized and Hapke-normalized data exist,
        // use intensity from the former and chromaticity from the latter.
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

void fillFace(const int order, const int pix, const uint8_t*const data,
              const ssize_t width, const ssize_t height, const ssize_t rowStride,
              const ssize_t channelsPerPixel, const QImage& lrocColorRef,
              const size_t refYMin, const size_t refYMax,
              const uint8_t*const smallData, double* outData)
{
    const unsigned nside = 1u << order;
    int ix, iy, face;
    healpix_nest2xyf(nside, pix, &ix, &iy, &face);

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
            fillPoint(longitude, latitude, pixelPosInOutData, data, width, height,
                      rowStride, lrocColorRef, refYMin, refYMax, smallData, outData);
        }
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
    QString lrocSRGBFileName;

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
        else if(arg == "--lroc-image")
        {
            GO_TO_PARAM();
            lrocSRGBFileName = argv[n];
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

    QString finalExt;
    if(imgFormat == "jpeg")
        finalExt = "jpg";
    else if(imgFormat == "webp" || imgFormat == "bmp" || imgFormat == "png" || imgFormat == "tiff")
        finalExt = imgFormat;
    else
        throw std::runtime_error("Unexpected output image format: " + imgFormat.toStdString());

    std::cerr << "Loading LROC-derived sRGB image... ";
    QImage lrocImg(lrocSRGBFileName);
    if(lrocImg.isNull())
        throw std::runtime_error(("Failed to load image \""+lrocSRGBFileName+"\"").toStdString());
    lrocImg = lrocImg.convertToFormat(QImage::Format_RGB888);
    if(lrocImg.isNull())
        throw std::runtime_error("Failed to convert the LROC-derived sRGB image to the working format");
    std::cerr << "done\n";

    std::cerr << "Loading input images...\n";
    constexpr ssize_t bigImgWidth = 27360*4;
    constexpr ssize_t bigImgHeight = bigImgWidth / 2;
    std::unique_ptr<uint8_t[]> bigData(new uint8_t[bigImgHeight*bigImgWidth]);

    const QString tiles[4][4] = {{"P900N0000-2250", "P900N0000-3150", "P900N0000-0450", "P900N0000-1350"},
                                 {"E300N2250"     , "E300N3150"     , "E300N0450"     , "E300N1350"     },
                                 {"E300S2250"     , "E300S3150"     , "E300S0450"     , "E300S1350"     },
                                 {"P900S0000-2250", "P900S0000-3150", "P900S0000-0450", "P900S0000-1350"}};
    ssize_t tileTopLineY = 0;
    for(const auto& row : tiles)
    {
        ssize_t tileLeftColX = 0;
        ssize_t tileHeight = -1;
        for(const auto& col : row)
        {
            const auto inFileName = QFileInfo(inDir + "/" + col + ".tiff").exists() ? col + ".tiff"
                                                                      : col + ".bmp";
            std::cerr << " " << inFileName.toStdString() << " ";
            const auto inFilePath = inDir + "/" + inFileName;
            QImage img(inFilePath);
            if(img.isNull())
                throw std::runtime_error(("Failed to open input image \""+inFilePath+"\"").toStdString());
            img = img.convertToFormat(QImage::Format_Grayscale8);
            if(img.isNull())
                throw std::runtime_error("Failed to convert input image to the working format");
            tileHeight = img.height();
            const auto tilePos = tileTopLineY*bigImgWidth + tileLeftColX;
            for(int j = 0; j < img.height(); ++j)
            {
                std::memcpy(bigData.get() + tilePos + j*bigImgWidth,
                            img.bits() + j*img.bytesPerLine(),
                            img.width());
            }
            tileLeftColX += img.width();
            std::cerr << "\n";
        }
        tileTopLineY += tileHeight;
    }
    std::cerr << "done\n";

    std::cerr << "Generating a lower-resolution version of the full map...";
    const ssize_t smallImgWidth  = bigImgWidth  / DOWNSAMP_FACTOR;
    const ssize_t smallImgHeight = bigImgHeight / DOWNSAMP_FACTOR;
    std::unique_ptr<uint8_t[]> smallData(new uint8_t[smallImgWidth*smallImgHeight]);
    for(ssize_t j = 0; j < smallImgHeight; ++j)
    {
        for(ssize_t i = 0; i < smallImgWidth; ++i)
        {
            double sum = 0;
            for(ssize_t x = 0; x < DOWNSAMP_FACTOR; ++x)
                for(ssize_t y = 0; y < DOWNSAMP_FACTOR; ++y)
                    sum += bigData[(j*DOWNSAMP_FACTOR+y)*bigImgWidth + (i*DOWNSAMP_FACTOR+x)];
            smallData[j*smallImgWidth+i] = sum / (DOWNSAMP_FACTOR*DOWNSAMP_FACTOR);
        }
    }
    std::cerr << "done\n";

    const int orderMax = std::ceil(std::log2(bigImgWidth / (4. * HIPS_TILE_SIZE * M_SQRT2)));

    if(!QDir().mkpath(outDir))
        throw std::runtime_error("Failed to create directory \""+outDir.toStdString()+'"');
    hipsSaveProperties(outDir, orderMax, imgFormat, surveyTitle, surveyType, description,
                       frame, obs_copyright, hips_copyright, creator, hipsStatus);

    // First create the tiles of the deepest level
    std::cerr << "Creating tiles of order " << orderMax << "...\n";
    {
        const int absolutePixMax = 12 * (1 << (2 * orderMax));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,orderMax,outDir,startTime,lrocImg,&smallData = std::as_const(smallData),
                     &bigData = std::as_const(bigData),colorChannel,
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
                fillFace(orderMax, pix, bigData.get(), bigImgWidth, bigImgHeight, bigImgWidth, channelsPerPixel,
                         lrocImg, refYmin, refYmax, smallData.get(), data.data());

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
                const auto fileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(hipsInitialExt);
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

    generateLowerOrderTiles(orderMax, outDir);
    convertTiles(finalExt, imgFormat, orderMax, outDir);
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
