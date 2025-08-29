#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <exception>
#include <tiffio.h>
#include <fitsio.h>
#include <Eigen/Dense>

using Float = double;
using Vec2d = Eigen::Matrix<Float, 2, 1>;
using Vec3d = Eigen::Matrix<Float, 3, 1>;
using Mat3d = Eigen::Matrix<Float, 3, 3>;
constexpr Float PI = Float(3.141592653589793238462643383279502884L);

constexpr Float sunR = 696000000;

class FITSError : public std::exception
{
    std::string msg;
public:
    FITSError(const std::string& prefix, const int status)
        : msg(prefix)
    {
        char error[FLEN_STATUS] = {};
        fits_get_errstatus(status, error);
        msg += ": ";
        msg += error;
    }
    const char* what() const noexcept override { return msg.c_str(); }
};

template<typename T> auto sqr(T x) { return x*x; }

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
Float sRGBTransferFunction(const Float c)
{
    const auto s = step(Float(0.0031308L),c);
    return s  * (Float(1.055L)*std::pow(c, 1/Float(2.4L))-Float(0.055L)) +
        (1-s) *  Float(12.92L)*c;
}

Mat3d Ry(const Float angle)
{
    Mat3d m;
    m = Eigen::AngleAxis<Float>(angle, -Vec3d::UnitY());
    return m;
}

Mat3d Rz(const Float angle)
{
    Mat3d m;
    m = Eigen::AngleAxis<Float>(angle, Vec3d::UnitZ());
    return m;
}

Vec2d obsXY(const Float sunLat, const Float sunLon,
            const Float carrLat, const Float carrLon,
            const Float distToObs)
{
    using namespace std;
    const Vec3d obsPos = distToObs * Vec3d(cos(carrLat)*cos(carrLon),
                                           cos(carrLat)*sin(carrLon),
                                           sin(carrLat));
    const Vec3d posRelToSun = sunR * Vec3d(cos(sunLat)*cos(sunLon),
                                           cos(sunLat)*sin(sunLon),
                                           sin(sunLat));
    const Vec3d posRelToObs = Rz(PI) * Ry(-carrLat) * Rz(-carrLon) * (posRelToSun - obsPos);

    if(posRelToSun.transpose() * obsPos < 0.) return {PI, PI/2};

    // The image is rotated by 180° from the north-at-top orientation,
    // so the signs here are opposite to the theoretical ones.
    return Vec2d(atan2(posRelToObs[1], posRelToObs[0]),
                 -asin(posRelToObs[2] / sqrt(posRelToObs.transpose() * posRelToObs)));
}

double fetchFromRect(double const*const data, const ssize_t stride, const ssize_t height,
                     const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = std::clamp(requestedX, ssize_t(0), stride-1);
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);
    return data[y*stride + x];
}

double sampleImg(const double*const data,
                 const ssize_t stride, const ssize_t height,
                 const Float x, const Float y)
{
    const auto floorX = std::floor(x);
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetchFromRect(data, stride, height, floorX  , floorY);
    const auto pTopRight    = fetchFromRect(data, stride, height, floorX+1, floorY);
    const auto pBottomLeft  = fetchFromRect(data, stride, height, floorX  , floorY+1);
    const auto pBottomRight = fetchFromRect(data, stride, height, floorX+1, floorY+1);
    // If at least one sample is invalid, discard the whole interpolant
    if(std::isnan(pTopLeft) || std::isnan(pTopRight) || std::isnan(pBottomLeft) || std::isnan(pBottomRight))
        return NAN;

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;
    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

double sampleImg(const double*const data, const int W, const int H, const int stride,
                 const double centerX, const double centerY,
                 const double scaleX, const double scaleY,
                 const Float sunLat, const Float sunLon,
                 const Float carrLat, const Float carrLon,
                 const Float sunObsDist, const Float maxAngularRadius)
{
    const Vec2d pos = obsXY(sunLat, sunLon, carrLat, carrLon, sunObsDist);
    if(pos.transpose() * pos > sqr(maxAngularRadius))
        return NAN;

    const auto i = centerX + pos[0]/scaleX;
    const auto j = H - 1 - (centerY + pos[1]/scaleY);
    if(0 <= i && i < W && 0 <= j && j < H)
        return sampleImg(data, stride, H, i, j);
    return NAN;
}

int usage(const char*const argv0, const int ret)
{
        std::cerr << "Usage: " << argv0 << 1+R"(
 input output
Options:
 --help                 Show this help message and quit
 --min NUM              Normalize output values so that NUM in the input becomes the smallest value of the output
 --max NUM              Normalize output values so that NUM in the input becomes the highest value of the output
 --crlt NUM             Override Carrington latitude
 --crln NUM             Override Carrington longitude
 --size NUM             Fraction of solar disk radius to use
 --uint8                Use 8 bpc in the output file (default)
 --uint16               Use 16 bpc in the output file
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    if(argc == 1)
        return usage(argv[0],1);

    bool useUInt16 = false;
    std::string infile;
    std::string outfile;
    int totalPositionalArgumentsFound = 0;
    double inputValueToOutputMax = 2;
    double inputValueToOutputMin = 0;
    double percentageOfSolarDiskToUse = 1;
    Float carrLon = NAN;
    Float carrLat = NAN;
    for(int n = 1; n < argc; ++n)
    {
#define REQUIRE_PARAM() do{                                                 \
            if(argc+1 == n)                                                 \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return usage(argv[0], 1);                                   \
            }}while(false)

        if(argv[n][0]!='-')
        {
            // Must be a positional argument
            switch(totalPositionalArgumentsFound)
            {
            case 0:
                infile = argv[n];
                break;
            case 1:
                outfile = argv[n];
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
        if(arg == "--help")
        {
            return usage(argv[0], 0);
        }
        else if(arg == "--min")
        {
            REQUIRE_PARAM();
            inputValueToOutputMin = std::stod(argv[++n]);
        }
        else if(arg == "--max")
        {
            REQUIRE_PARAM();
            inputValueToOutputMax = std::stod(argv[++n]);
        }
        else if(arg == "--crlt")
        {
            REQUIRE_PARAM();
            carrLat = std::stod(argv[++n]);
        }
        else if(arg == "--crln")
        {
            REQUIRE_PARAM();
            carrLon = std::stod(argv[++n]);
        }
        else if(arg == "--size")
        {
            REQUIRE_PARAM();
            percentageOfSolarDiskToUse = std::stod(argv[++n]);
        }
        else if(arg == "--uint16")
        {
            useUInt16 = true;
        }
        else if(arg == "--uint8")
        {
            useUInt16 = false;
        }
        else
        {
            std::cerr << "Unknown option " << arg << "\n";
            return usage(argv[0], 1);
        }
    }
#undef REQUIRE_PARAM

    fitsfile* fits;
    int status = 0;
    if(fits_open_diskfile(&fits, infile.c_str(), READONLY, &status))
        throw FITSError("Failed to open \""+infile+'"', status);
    int numHDUs;
    fits_get_num_hdus(fits, &numHDUs, &status);
    if(status) throw FITSError("Failed to get HDU count in the FITS file", status);
    if(numHDUs <= 0) throw std::runtime_error("No images in the FITS file");
    int bitpix, naxis;
    long naxes[3];
    int width = 0, height = 0;
    std::unique_ptr<double[]> data;
    for(int hdu = 1; hdu <= numHDUs; ++hdu)
    {
        fits_movabs_hdu(fits, hdu, nullptr, &status);
        if(status) throw FITSError("Failed to move to HDU #"+std::to_string(hdu)+"", status);
        fits_get_img_param(fits, std::size(naxes), &bitpix, &naxis, naxes, &status);
        if(status) throw FITSError("Failed to get image parameters of HDU #"+std::to_string(hdu)+"", status);
        if(naxis != 2) continue;
        width  = naxes[0];
        height = naxes[1];
        if(bitpix != SHORT_IMG)
            throw std::runtime_error("FITS image type isn't SHORT, instead bitpix="+std::to_string(bitpix)+", this is not supported");
        data.reset(new double[width*height]);
    }
    if(!data) throw std::runtime_error("Failed to find a usable image in the FITS file");
    const int stride = width;
    long fpixel[2] = {1,1};
    double nullVal = NAN;
    // Read the image from top to bottom, taking into account that cfitsio's fpixel[1]==1 corrsponds to the bottom scanline.
    for(int j = 0; j < height; ++j)
    {
        fpixel[1] = height - j;
        fits_read_pix(fits, TDOUBLE, fpixel, width, &nullVal, data.get() + width * j, nullptr, &status);
        if(status) throw FITSError("Failed to read pixels from the FITS image", status);
    }

    Float sunObsDist = NAN;
    Float centerX = NAN, centerY = NAN;
    Float scaleX = NAN, scaleY = NAN;

    if(std::isnan(carrLon))
        fits_read_key(fits, TDOUBLE, "CRLN_OBS", &carrLon, nullptr, &status); if(status) carrLon = NAN; status = 0;
    if(std::isnan(carrLat))
    fits_read_key(fits, TDOUBLE, "CRLT_OBS", &carrLat, nullptr, &status); if(status) carrLat = NAN; status = 0;
    fits_read_key(fits, TDOUBLE, "DSUN_OBS", &sunObsDist, nullptr, &status); if(status) sunObsDist = NAN; status = 0;
    fits_read_key(fits, TDOUBLE, "CRPIX1", &centerX, nullptr, &status); if(status) centerX = NAN; status = 0;
    fits_read_key(fits, TDOUBLE, "CRPIX2", &centerY, nullptr, &status); if(status) centerY = NAN; status = 0;
    fits_read_key(fits, TDOUBLE, "CDELT1", &scaleX, nullptr, &status); if(status) scaleX = NAN; status = 0;
    fits_read_key(fits, TDOUBLE, "CDELT2", &scaleY, nullptr, &status); if(status) scaleY = NAN; status = 0;

    if(std::isnan(carrLon) || std::isnan(carrLat) || std::isnan(sunObsDist) ||
       std::isnan(centerX) || std::isnan(centerY) || std::isnan(scaleX) || std::isnan(scaleY))
    {
        std::cerr << "Some keywords aren't defined in the FITS header, trying to read them from the sidecar file...\n";
        std::string sidecar;
        if(infile.size() >= 5 && infile[infile.size() - 5] == '.' && infile.substr(infile.size() - 5) == ".fits")
            sidecar = infile.substr(0, infile.size() - 5) + ".txt";
        else
            throw std::runtime_error("Input file name is strange (should end in .fits), can't find out sidecar name");
        std::ifstream s(sidecar);
        if(!s) throw std::runtime_error("Failed to open sidecar file "+sidecar);
        s.exceptions(std::ios_base::badbit);
        while(s)
        {
            std::string name;
            s >> name;
            if(name == "DSUN_OBS" && std::isnan(sunObsDist))
                s >> sunObsDist;
            else if(name == "CRLN_OBS" && std::isnan(carrLon))
                s >> carrLon;
            else if(name == "CRLT_OBS" && std::isnan(carrLat))
                s >> carrLat;
            else if(name == "CRPIX1" && std::isnan(centerX))
                s >> centerX;
            else if(name == "CRPIX2" && std::isnan(centerY))
                s >> centerY;
            else if(name == "CDELT1" && std::isnan(scaleX))
                s >> scaleX;
            else if(name == "CDELT2" && std::isnan(scaleY))
                s >> scaleY;
        }
    }
    if(std::isnan(carrLon)) throw std::runtime_error("CRLN_OBS field is missing in the sidecar file");
    if(std::isnan(carrLat)) throw std::runtime_error("CRLT_OBS field is missing in the sidecar file");
    if(std::isnan(sunObsDist)) throw std::runtime_error("DSUN_OBS field is missing in the sidecar file");
    if(std::isnan(centerX)) throw std::runtime_error("CRPIX1 field is missing in the sidecar file");
    if(std::isnan(centerY)) throw std::runtime_error("CRPIX2 field is missing in the sidecar file");
    if(std::isnan(scaleX)) throw std::runtime_error("CDELT1 field is missing in the sidecar file");
    if(std::isnan(scaleY)) throw std::runtime_error("CDELT2 field is missing in the sidecar file");
    carrLon *= PI / 180;
    carrLat *= PI / 180;
    scaleX *= PI / 180 / 3600;
    scaleY *= PI / 180 / 3600;
    centerX -= 1;
    centerY -= 1;
    std::cerr << "All keywords read successfully\n";

    double minVal = INFINITY, maxVal = -INFINITY;
    const ssize_t N = ssize_t(width) * height;
    for(ssize_t i = 0; i < N; ++i)
    {
        const auto v = data[i];
        if(std::isnan(v)) continue;
        if(v > maxVal) maxVal = v;
        if(v < minVal) minVal = v;
    }
    std::cerr << "Minimum value: " << minVal << ", maximum value: " << maxVal << "\n";

    const Float sunAngR = asin(sunR/sunObsDist);
    const Float maxAngularRadius = sunAngR * percentageOfSolarDiskToUse;

    constexpr ssize_t outH = 4096, outW = 2*outH, outSPP = 3;

    TIFF*const tif = TIFFOpen(outfile.c_str(), "w");
    if(!tif) throw std::runtime_error("Failed to open file \"" + outfile + '"');
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, outW);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, outH);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, outSPP);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    auto process = [&](const auto& type)
    {
        using OutType = std::remove_cv_t<std::remove_reference_t<decltype(type)>>;
        std::unique_ptr<OutType[]> outScanLine(new OutType[outW * outSPP]);
        constexpr ssize_t outBPS = sizeof outScanLine[0] * 8;
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, outBPS);

        constexpr auto OUT_VAL_MAX = std::numeric_limits<std::remove_reference_t<decltype(outScanLine[0])>>::max();
        for(ssize_t j = 0; j < outH; ++j)
        {
            const Float sunLat = -(Float(j) - outH / Float(2)) / outH * PI;
            for(ssize_t i = 0; i < outW; ++i)
            {
                const Float sunLon = (Float(i + 0.5) - outW / Float(2)) / outW * (2 * PI);
                const auto samp = sampleImg(data.get(), width, height, stride, centerX, centerY,
                                            scaleX, scaleY, sunLat, sunLon, carrLat, carrLon,
                                            sunObsDist, maxAngularRadius);
                const ssize_t index = i * outSPP;
                const auto val = std::isnan(samp) ? 0 :
                    (samp - inputValueToOutputMin) / (inputValueToOutputMax - inputValueToOutputMin);
                const auto valR = val * 1;
                const auto valG = val * 0.890095;
                const auto valB = val * 0.859308;
                outScanLine[index + 0] = OUT_VAL_MAX * sRGBTransferFunction(std::clamp(valR, 0., 1.));
                outScanLine[index + 1] = OUT_VAL_MAX * sRGBTransferFunction(std::clamp(valG, 0., 1.));
                outScanLine[index + 2] = OUT_VAL_MAX * sRGBTransferFunction(std::clamp(valB, 0., 1.));
            }
            if(TIFFWriteScanline(tif, outScanLine.get(), j, 0) != 1)
                throw std::runtime_error("Failed to write scan line "+std::to_string(j + 1));
        }
    };
    if(useUInt16)
        process(uint16_t(0));
    else
        process(uint8_t(0));

    TIFFClose(tif);
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
