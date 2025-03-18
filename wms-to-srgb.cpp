#include <iostream>
#include <QImage>

namespace
{

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    const auto s = step(0.0031308,c);
    return s  * (1.055*std::pow(c, 1/2.4)-0.055) +
        (1-s) *  12.92*c;
}

double wmsToRed(const double x)
{
    return 3.66897675367758 + x*(5.17181626686068 + x*(37.4567310735643 + x*(-66.1407972224146 + 36.9574791372878*x)));
}
double wmsToGreen(const double x)
{
    return 2.72459214650356 + x*(4.07386797773892 + x*(26.3036294920035 + x*(-44.8836716108057 + 25.0812660136267*x)));
}
double wmsToBlue(const double x)
{
    return 2.09428707617804 + x*(2.89641075271338 + x*(18.7590518760633 + x*(-29.7670405560604 + 16.4039924091487*x)));
}

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    if(!ret)
    {
        s << "This program tries to convert the colors of the LROC WAC Normalized Reflectance layer "
             "from https://wms.lroc.asu.edu/lroc to sRGB, aiming at getting much the same result as "
             "lroc-raw-to-srgb does with the actual per-wavelength data.\n";
    }
    s << "Usage: " << argv0 << "{options...} inputFile outputFile";
    s << R"(
Options:
 -h, --help                 This help message plus a little description
 --value-scale N            Scale albedo by N before saving as color
)";
    return ret;
}
}

int main(int argc, char** argv)
try
{
    QString inFileName;
    QString outFileName;
    constexpr double defaultScale = 105;
    double valueScale = 1 / defaultScale;

    int totalPositionalArgumentsFound = 0;
    for(int n = 1; n < argc; ++n)
    {
        if(argv[n][0]!='-')
        {
            // Must be a positional argument
            switch(totalPositionalArgumentsFound)
            {
            case 0:
                inFileName = argv[n];
                break;
            case 1:
                outFileName = argv[n];
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

        if(arg == "--value-scale")
        {
            GO_TO_PARAM();
            valueScale = std::stod(argv[n]) / defaultScale;
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(inFileName.isEmpty())
    {
        std::cerr << "Input filename not specified\n";
        return usage(argv[0], 1);
    }
    if(outFileName.isEmpty())
    {
        std::cerr << "Output file name not specified\n";
        return usage(argv[0], 1);
    }

    std::cerr << "Loading image...\n";
    QImage img(inFileName);
    if(img.isNull())
    {
        std::cerr << "Failed to load input image\n";
        return 1;
    }
    img = img.convertToFormat(QImage::Format_RGB888);
    const uchar*const data = img.bits();
    constexpr int bytesPerColor = 3;
    const ssize_t rowStride = img.bytesPerLine() / bytesPerColor;
    const ssize_t width = img.width();
    const ssize_t height = img.height();

    std::cerr << "Converting data to colors...\n";
    std::vector<double> out(rowStride*height*bytesPerColor);
    for(size_t n = 0; n < out.size(); n += bytesPerColor)
    {
        const double v = data[n] / 255.;
        out[n+0] = wmsToRed(v);
        out[n+1] = wmsToGreen(v);
        out[n+2] = wmsToBlue(v);
    }

    std::vector<uint8_t> outImgData(rowStride*height*bytesPerColor);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(std::abs(out[n]) * valueScale, 0., 1.));
    const QImage outImg(outImgData.data(), width, height, rowStride*bytesPerColor, QImage::Format_RGB888);
    std::cerr << "Saving output file with dimensions " << width << " × " << height << "... ";
    if(!outImg.save(outFileName))
    {
        std::cerr << "Failed to save output file " << outFileName.toStdString() << "\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
