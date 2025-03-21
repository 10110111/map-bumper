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

const float* getSurroundingColor(const float*const data,
                                  const ssize_t width, const ssize_t height, const int bytesPerColor,
                                  const ssize_t i0, const ssize_t j0)
{
    struct Color
    {
        float r,g,b;
        Color() : r(0), g(0), b(0) {}
        Color(const float* c) : r(c[0]), g(c[1]), b(c[2]) {}
    };

    static std::vector<Color> colors;
    colors.clear();
    static std::vector<ssize_t> radii;
    radii.clear();
    static constexpr int directions[][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1},{-1,1},{1,-1},
                                            {2,1},{1,2},{-1,2},{-2,1},{-2,-1},{-1,-2},{1,-2},{2,-1}};
    for(const auto dir : directions)
    {
        const double dirNorm = std::sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
        ssize_t radius = 0;
        ssize_t i, j;
        do
        {
            ++radius;
            i = i0+radius*dir[0];
            j = j0+radius*dir[1];
            if(!(0 <= i && i < width && 0 <= j && j < height))
                goto unusableDirection;
            if(data[(j*width+i)*bytesPerColor+1] > 0)
                break; // found a valid value
        } while(true);
        colors.emplace_back(data + ((j0+radius*dir[1])*width+(i0+radius*dir[0]))*bytesPerColor);
        if(colors.back().g <= 0) radius = 0;
        radii.push_back(radius * dirNorm);

unusableDirection:
        continue;
    }

    static Color out;
    out = {};
    double commonCoef = 0;
    for(const auto r : radii)
        commonCoef += 1. / r;
    for(unsigned n = 0; n < colors.size(); ++n)
    {
        out.r += colors[n].r / radii[n];
        out.g += colors[n].g / radii[n];
        out.b += colors[n].b / radii[n];
    }
    out.r /= commonCoef;
    out.g /= commonCoef;
    out.b /= commonCoef;
    return &out.r;
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
    setenv("QT_IMAGEIO_MAXALLOC","2048", false);

    QString inFileName;
    QString outFileName;
    constexpr double defaultScale = 105;
    double valueScale = 1 / defaultScale;
    double badLevel = 1. / 255.;
    enum class MarkBadMode
    {
        Magenta,
        Surroundings,
        None
    };
    MarkBadMode markBadMode = MarkBadMode::None;

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
        else if(arg == "--bad-level")
        {
            GO_TO_PARAM();
            badLevel = std::stod(argv[n]) / 255.;
        }
        else if(arg == "--fill-bad")
        {
            markBadMode = MarkBadMode::Surroundings;
        }
        else if(arg == "--bad-to-magenta")
        {
            markBadMode = MarkBadMode::Magenta;
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
    std::vector<float> out(rowStride*height*bytesPerColor);
    for(size_t n = 0; n < out.size(); n += bytesPerColor)
    {
        const auto j = n / img.bytesPerLine();
        const double v = data[n] / 255.;
        const bool good = v > badLevel || (1520+1300 < j && j < 1520+9570);
        if(good || markBadMode == MarkBadMode::None)
        {
            out[n+0] = wmsToRed(v);
            out[n+1] = wmsToGreen(v);
            out[n+2] = wmsToBlue(v);
        }
        else
        {
            // Mark with magenta
            out[n+0] = 1 / valueScale;
            out[n+1] = 0;
            out[n+2] = 1 / valueScale;
        }
    }

    if(markBadMode == MarkBadMode::Surroundings)
    {
        std::cerr << "Replacing bad pixels with some surrounding colors...\n";
        for(ssize_t j = 0; j < height; ++j)
        {
            const ssize_t linePos = j*width;
            for(ssize_t i = 0; i < width; ++i)
            {
                const auto pos = (linePos + i)*bytesPerColor;
                if(out[pos+1] == 0)
                {
                    const auto*const surroundingColor = getSurroundingColor(out.data(), width, height, bytesPerColor, i, j);
                    assert(surroundingColor);
                    if(!surroundingColor) continue;
                    // Mark fixed colors with a negative sign to prevent using them as surroundings.
                    // The negativity will be removed by abs() when saving the image.
                    out[pos+0] = -surroundingColor[0];
                    out[pos+1] = -surroundingColor[1];
                    out[pos+2] = -surroundingColor[2];
                }
            }
        }
        // Blur the replaced pixels a bit to better fit into the surroundings
        for(ssize_t j = 0; j < height; ++j)
        {
            const ssize_t linePos = j*width;
            for(ssize_t i = 0; i < width; ++i)
            {
                const auto pos = (linePos + i)*bytesPerColor;
                if(out[pos+1] < 0)
                {
                    using std::abs;
                    for(int c = 0; c < 3; ++c)
                    {
                        const auto left = std::max(i-1, ssize_t(0));
                        const auto right = std::min(i+1, width-1);
                        const auto top = std::max(j-1, ssize_t(0));
                        const auto bottom = std::min(j+1, height-1);
                        const auto B = bytesPerColor;
                        out[pos+c] = abs(out[(top*width+left)*B+c]) + abs(out[(top*width+i)*B+c]) + abs(out[(top*width+right)*B+c])+
                                     abs(out[((j+0)*width+left)*B+c]) + abs(out[((j+0)*width+i)*B+c]) + abs(out[((j+0)*width+right)*B+c])+
                                     abs(out[(bottom*width+left)*B+c]) + abs(out[(bottom*width+i)*B+c]) + abs(out[(bottom*width+right)*B+c]);
                        out[pos+c] /= 9;
                    }
                }
            }
        }
    }

    std::vector<uint8_t> outImgData(rowStride*height*bytesPerColor);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(std::abs(out[n]) * valueScale, 0., 1.));

    // Clear some RAM for the saving procedure
    out.clear(); out.shrink_to_fit();

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
