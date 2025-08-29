#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <exception>
#include <tiffio.h>
#include <fitsio.h>

using Float = double;
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

    constexpr ssize_t outSPP = 1;
    using OutType = uint16_t;
    std::unique_ptr<OutType[]> outScanLine(new OutType[width * outSPP]);
    constexpr ssize_t outBPS = sizeof outScanLine[0] * 8;

    TIFF*const tif = TIFFOpen(outfile.c_str(), "w");
    if(!tif) throw std::runtime_error("Failed to open file \"" + outfile + '"');
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, outBPS);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, outSPP);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    constexpr auto OUT_VAL_MAX = std::numeric_limits<std::remove_reference_t<decltype(outScanLine[0])>>::max();
    for(ssize_t j = 0; j < height; ++j)
    {
        for(ssize_t i = 0; i < width; ++i)
        {
            const auto samp = data[stride * j + i];
            const ssize_t index = i * outSPP;
            const auto val = std::isnan(samp) ? 0 : samp / inputValueToOutputMax;
            outScanLine[index] = OUT_VAL_MAX * sRGBTransferFunction(std::clamp(val, 0., 1.));
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
