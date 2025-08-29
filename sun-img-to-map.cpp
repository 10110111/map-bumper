#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <tiffio.h>
#include <Eigen/Dense>

using Vec2d = Eigen::Vector2d;
using Vec3d = Eigen::Vector3d;
using Mat3d = Eigen::Matrix3d;

constexpr int SAMP_PER_PX = 2;
constexpr double sunR = 696000000.;

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    const auto s = step(0.0031308,c);
    return s  * (1.055*std::pow(c, 1/2.4)-0.055) +
        (1-s) *  12.92*c;
}

Mat3d Ry(const double angle)
{
    Mat3d m;
    m = Eigen::AngleAxisd(angle, -Vec3d::UnitY());
    return m;
}

Mat3d Rz(const double angle)
{
    Mat3d m;
    m = Eigen::AngleAxisd(angle, Vec3d::UnitZ());
    return m;
}

Vec2d obsXY(const double sunLat, const double sunLon,
            const double carrLat, const double carrLon,
            const double distToObs)
{
    using namespace std;
    const Vec3d obsPos = distToObs * Vec3d(cos(carrLat)*cos(carrLon),
                                           cos(carrLat)*sin(carrLon),
                                           sin(carrLat));
    const Vec3d posRelToSun = sunR * Vec3d(cos(sunLat)*cos(sunLon),
                                           cos(sunLat)*sin(sunLon),
                                           sin(sunLat));
    const Vec3d posRelToObs = Rz(M_PI) * Ry(-carrLat) * Rz(-carrLon) * (posRelToSun - obsPos);

    if(posRelToSun.transpose() * obsPos < 0.) return {M_PI, M_PI/2};

    return Vec2d(atan2(posRelToObs[1], posRelToObs[0]),
                 asin(posRelToObs[2] / sqrt(posRelToObs.transpose() * posRelToObs)));
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

Vec2d sampleImg(const uint16_t*const data, const int W, const int H, const int stride,
                const int imgWidth, const int imgHeight,
                const int marginLeft, const int marginRight, const int marginTop, const int marginBottom,
                const double sunLat, const double sunLon,
                const double carrLat, const double carrLon,
                const double sunObsDist, const double sunAngR)
{
    const Vec2d pos = obsXY(sunLat, sunLon, carrLat, carrLon, sunObsDist);
    // FIXME: interpolate linearly instead of to nearest neighbour
    const ssize_t i = std::lround(marginLeft + imgWidth /2. + pos[0]/sunAngR*imgWidth/2.);
    const ssize_t j = std::lround(marginTop  + imgHeight/2. + pos[1]/sunAngR*imgHeight/2.);
    if(marginLeft <= i && i < W - marginRight &&
       marginTop  <= j && j < H - marginBottom)
    {
        return {data[j * stride + i * SAMP_PER_PX + 0],
                data[j * stride + i * SAMP_PER_PX + 1]};
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

    // TODO: read these from a sidecar .txt file
    const double carrLon = 207.310318 * M_PI / 180;
    const double carrLat = 7.007775 * M_PI / 180;
    const double sunObsDist = 151252300213.76;

    const double sunAngR = asin(sunR/sunObsDist);

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
        const double sunLat = (double(j) - outH / 2.) / outH * M_PI;
        for(ssize_t i = 0; i < outW; ++i)
        {
            const double sunLon = (double(i) - outW / 2.) / outW * (2 * M_PI);
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
