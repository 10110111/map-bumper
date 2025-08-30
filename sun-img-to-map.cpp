#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <tiffio.h>
#include <Eigen/Dense>

using Float = double;
using Vec2d = Eigen::Matrix<Float, 2, 1>;
using Vec3d = Eigen::Matrix<Float, 3, 1>;
using Mat3d = Eigen::Matrix<Float, 3, 3>;
constexpr Float PI = Float(3.141592653589793238462643383279502884L);

constexpr int SAMP_PER_PX = 2;
constexpr Float sunR = 696000000;

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

void getMargins(const uint16_t*const data, const int W, const int H,
                const int stride,
                int& left, int& right, int& top, int& bottom)
{
    left = 0;
    for(int n = 0; n < W/2; ++n)
        if(data[H/2 * stride + n * SAMP_PER_PX + 1] == 0)
            ++left;
    right = 0;
    for(int n = W - 1; n >= W/2; --n)
        if(data[H/2 * stride + n * SAMP_PER_PX + 1] == 0)
            ++right;
    top = 0;
    for(int n = 0; n < H/2; ++n)
        if(data[n * stride + W/2 * SAMP_PER_PX + 1] == 0)
            ++top;
    bottom = 0;
    for(int n = H - 1; n >= H/2; --n)
        if(data[n * stride + W/2 * SAMP_PER_PX + 1] == 0)
            ++bottom;
}

Vec2d fetchFromRect(uint16_t const*const data, const ssize_t stride, const ssize_t height,
                    const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = std::clamp(requestedX, ssize_t(0), stride-1);
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);
    return {data[y*stride + x * SAMP_PER_PX + 0],
            data[y*stride + x * SAMP_PER_PX + 1]};
}

Vec2d sampleImg(const uint16_t*const data,
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
    if(!pTopLeft[1] || !pTopRight[1] || !pBottomLeft[1] || !pBottomRight[1])
        return {NAN, 0};

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;
    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

Vec2d sampleImg(const uint16_t*const data, const int W, const int H, const int stride,
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
    return {NAN, 0};
}

int main(int argc, char** argv)
try
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " input output\n";
        return 1;
    }

    const std::string infile=argv[1];
    const std::string outfile=argv[2];
    if(infile.size() < 4)
        throw std::runtime_error("Input file name is strange (should end in .tiff or .tif), can't find out sidecar name");
    std::string sidecar;
    if(infile[infile.size() - 4] == '.' && infile.substr(infile.size() - 4) == ".tif")
        sidecar = infile.substr(0, infile.size() - 4) + ".txt";
    else if(infile.size() >= 5 && infile[infile.size() - 5] == '.' && infile.substr(infile.size() - 5) == ".tiff")
        sidecar = infile.substr(0, infile.size() - 5) + ".txt";
    else
        throw std::runtime_error("Input file name is strange (should end in .tiff or .tif), can't find out sidecar name");
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

    TIFF* tif = TIFFOpen(infile.c_str(), "r");
    if(!tif)
        throw std::runtime_error("Failed to open file \""+infile+'"');
    uint32_t width = 0;
    if(!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width))
        throw std::runtime_error("Failed to get image width");
    uint32_t height = 0;
    if(!TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height))
        throw std::runtime_error("Failed to get image height");
    uint32_t bpp = 0;
    if(!TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bpp))
        throw std::runtime_error("Failed to get image height");
    uint32_t planarConf = 0;
    if(!TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &planarConf))
        throw std::runtime_error("Failed to get image height");
    if(planarConf != PLANARCONFIG_CONTIG)
        throw std::runtime_error("TIFF planar config is not contiguous, this is not supported");
    uint32_t sampFmt = 0;
    if(!TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sampFmt))
        throw std::runtime_error("Failed to get sample format");
    if(sampFmt != SAMPLEFORMAT_UINT)
        throw std::runtime_error("Sample format is not uint: it's "+std::to_string(sampFmt)+", it's not supported");
    uint32_t spp = 0;
    if(!TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp))
        throw std::runtime_error("Failed to get number of samples per pixel");
    if(spp != 2)
        throw std::runtime_error("The image has "+std::to_string(spp)+" samples per pixel, it's not supported");
    uint32_t phot = 0;
    if(!TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &phot))
        throw std::runtime_error("Failed to get photometric interpretation");
    if(phot != PHOTOMETRIC_MINISBLACK)
        throw std::runtime_error("Unsupported photometric interpretation "+std::to_string(phot));

    const auto stride = TIFFScanlineSize(tif);
    const std::unique_ptr<uint16_t[]> data(new uint16_t[stride*height]);
    for(unsigned row = 0; row < height; ++row)
    {
        if(!TIFFReadScanline(tif, &data[row * stride], row))
            throw std::runtime_error("Failed to read scanline "+std::to_string(row+1));
    }
    TIFFClose(tif);

    int marginLeft, marginRight, marginTop, marginBottom;
    getMargins(data.get(), width, height, stride, marginLeft, marginRight, marginTop, marginBottom);
    const int imgWidth = width - marginLeft - marginRight;
    const int imgHeight = height - marginTop - marginBottom;

    const Float sunAngR = asin(sunR/sunObsDist);

    constexpr ssize_t outH = 4096, outW = 2*outH, outBytesPP = 1, outSPP = 1;
    std::unique_ptr<uint8_t[]> outScanLine(new uint8_t[outW * outBytesPP]);

    tif = TIFFOpen(outfile.c_str(), "w");
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
            outScanLine[index] = samp[1] ? samp[0] / 65535. * 255. : 0;
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
