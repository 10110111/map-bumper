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

void getMargins(const double*const data, const int W, const int H,
                const int stride,
                int& left, int& right, int& top, int& bottom)
{
    left = 0;
    for(int n = 0; n < W/2; ++n)
        if(std::isnan(data[H/2 * stride + n]))
            ++left;
    right = 0;
    for(int n = W - 1; n >= W/2; --n)
        if(std::isnan(data[H/2 * stride + n]))
            ++right;
    top = 0;
    for(int n = 0; n < H/2; ++n)
        if(std::isnan(data[n * stride + W/2]))
            ++top;
    bottom = 0;
    for(int n = H - 1; n >= H/2; --n)
        if(std::isnan(data[n * stride + W/2]))
            ++bottom;
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
                 const int imgWidth, const int imgHeight,
                 const int marginLeft, const int marginRight, const int marginTop, const int marginBottom,
                 const Float sunLat, const Float sunLon,
                 const Float carrLat, const Float carrLon,
                 const Float sunObsDist, const Float sunAngR)
{
    const Vec2d pos = obsXY(sunLat, sunLon, carrLat, carrLon, sunObsDist);
    const auto i = marginLeft + imgWidth /2. + pos[0]/sunAngR*imgWidth/2.;
    const auto j = marginTop  + imgHeight/2. - pos[1]/sunAngR*imgHeight/2.;
    if(marginLeft <= i && i < W - marginRight &&
       marginTop  <= j && j < H - marginBottom)
    {
        return sampleImg(data, stride, H, i, j);
    }
    return NAN;
}

int usage(const char*const argv0, const int ret)
{
        std::cerr << "Usage: " << argv0 << 1+R"(
 input output
Options:
 --help                 Show this help message and quit
 -m,--max NUM           Normalize output values so that NUM in the input becomes the highest value of the output
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    if(argc == 1)
        return usage(argv[0],1);

    std::string infile;
    std::string outfile;
    int totalPositionalArgumentsFound = 0;
    double inputValueToOutputMax = 2;
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
        else if(arg == "-m" || arg == "--max")
        {
            REQUIRE_PARAM();
            inputValueToOutputMax = std::stod(argv[++n]);
        }
        else
        {
            std::cerr << "Unknown option " << arg << "\n";
            return usage(argv[0], 1);
        }
    }
#undef REQUIRE_PARAM

    std::string sidecar;
    if(infile.size() >= 5 && infile[infile.size() - 5] == '.' && infile.substr(infile.size() - 5) == ".fits")
        sidecar = infile.substr(0, infile.size() - 5) + ".txt";
    else
        throw std::runtime_error("Input file name is strange (should end in .fits), can't find out sidecar name");
    Float carrLon = NAN;
    Float carrLat = NAN;
    Float sunObsDist = NAN;
    {
        std::ifstream s(sidecar);
        if(!s) throw std::runtime_error("Failed to open sidecar file "+sidecar);
        s.exceptions(std::ios_base::badbit);
        while(s)
        {
            std::string name;
            s >> name;
            if(name == "DSUN_OBS")
                s >> sunObsDist;
            else if(name == "CRLN_OBS")
                s >> carrLon;
            else if(name == "CRLT_OBS")
                s >> carrLat;
        }
        if(std::isnan(carrLon)) throw std::runtime_error("CRLN_OBS field is missing in the sidecar file");
        if(std::isnan(carrLat)) throw std::runtime_error("CRLT_OBS field is missing in the sidecar file");
        if(std::isnan(sunObsDist)) throw std::runtime_error("DSUN_OBS field is missing in the sidecar file");
        carrLon *= PI / 180;
        carrLat *= PI / 180;
    }

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

    int marginLeft, marginRight, marginTop, marginBottom;
    getMargins(data.get(), width, height, stride, marginLeft, marginRight, marginTop, marginBottom);
    const int imgWidth = width - marginLeft - marginRight;
    const int imgHeight = height - marginTop - marginBottom;

    const Float sunAngR = asin(sunR/sunObsDist);

    constexpr ssize_t outH = 4096, outW = 2*outH, outBytesPP = 1, outSPP = 1;
    std::unique_ptr<uint8_t[]> outScanLine(new uint8_t[outW * outBytesPP]);

    TIFF*const tif = TIFFOpen(outfile.c_str(), "w");
    if(!tif) throw std::runtime_error("Failed to open file \"" + outfile + '"');
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, outW);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, outH);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, outBytesPP * 8 / outSPP);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, outSPP);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    for(ssize_t j = 0; j < outH; ++j)
    {
        const Float sunLat = -(Float(j) - outH / Float(2)) / outH * PI;
        for(ssize_t i = 0; i < outW; ++i)
        {
            const Float sunLon = (Float(i) - outW / Float(2)) / outW * (2 * PI);
            const auto samp = sampleImg(data.get(), width, height, stride, imgWidth, imgHeight,
                                        marginLeft, marginRight, marginTop, marginBottom,
                                        sunLat, sunLon, carrLat, carrLon, sunObsDist, sunAngR);
            const ssize_t index = i * outBytesPP;
            const auto val = std::isnan(samp) ? 0 : samp / inputValueToOutputMax;
            outScanLine[index] =  255 * sRGBTransferFunction(std::clamp(val, 0., 1.));
        }
        if(TIFFWriteScanline(tif, outScanLine.get(), j, 0) != 1)
            throw std::runtime_error("Failed to write scan line "+std::to_string(j + 1));
    }

    TIFFClose(tif);
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
