#include <iostream>

#include <QFile>
#include <QImage>
#include <QByteArray>
#include <QRegularExpression>

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    const auto s = step(0.0031308,c);
    return s  * (1.055*std::pow(c, 1/2.4)-0.055) +
        (1-s) *  12.92*c;
}

const double* getSurroundingColor(const double*const data,
                                  const ssize_t width, const ssize_t height, const int bytesPerColor,
                                  const ssize_t i0, const ssize_t j0)
{
    struct Color
    {
        double r,g,b;
        Color(const double* c) : r(c[0]), g(c[1]), b(c[2]) {}
    };
    static std::vector<Color> colors;

    const auto maxRadius = std::max(width, height);
    for(int radius = 1; radius < maxRadius; ++radius)
    {
        colors.clear();
        {
            ssize_t j = std::clamp(j0-radius, ssize_t(0), height-1);
            const ssize_t iMin = std::clamp(i0-radius, ssize_t(0), width-1);
            const ssize_t iMax = std::clamp(i0+radius, ssize_t(0), width-1);
            for(ssize_t i = iMin; i <= iMax; ++i)
            {
                const auto pos = (j*width + i)*bytesPerColor;
                if(data[pos+1] > 0) colors.push_back(data+pos);
            }
            j = std::clamp(j0+radius, ssize_t(0), height-1);
            for(ssize_t i = iMin; i <= iMax; ++i)
            {
                const auto pos = (j*width + i)*bytesPerColor;
                if(data[pos+1] > 0) colors.push_back(data+pos);
            }
        }
        {
            ssize_t i = std::clamp(i0-radius, ssize_t(0), width-1);
            const ssize_t jMin = std::clamp(j0-radius, ssize_t(0), height-1);
            const ssize_t jMax = std::clamp(j0+radius, ssize_t(0), height-1);
            for(ssize_t j = jMin; j <= jMax; ++j)
            {
                const auto pos = (j*width + i)*bytesPerColor;
                if(data[pos+1] > 0) colors.push_back(data+pos);
            }
            i = std::clamp(i0+radius, ssize_t(0), width-1);
            for(ssize_t j = jMin; j <= jMax; ++j)
            {
                const auto pos = (j*width + i)*bytesPerColor;
                if(data[pos+1] > 0) colors.push_back(data+pos);
            }
        }
        if(colors.size() < 4)
        {
            // Too few valid colors may be undetected bads
            colors.clear();
            continue;
        }

        // Return the median of a simplified brightness
        std::sort(colors.begin(), colors.end(), [](const auto& a, const auto& b) { return a.r+a.g+a.b < b.r+b.g+b.b; });
        return &colors[std::min(colors.size() / 2, colors.size())].r;
    }
    return nullptr;
}

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} inDir sector outFile";
    s << R"(
Here sector is a part of the input filenames denoting map sector, e.g. E350N0450.
Options:
 -h, --help                 This help message
 --value-scale N            Scale albedo by N before saving as color
 --bad-to-magenta           Fill bad pixels with magenta
 --fill-bad                 Fill bad pixels with surrounding colors
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    std::string inDir;
    std::string sector;
    std::string outFileName;
    constexpr double defaultScale = 105;
    double valueScale = 1 / defaultScale;
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
                inDir = argv[n];
                break;
            case 1:
                sector = argv[n];
                break;
            case 2:
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
    if(inDir.empty())
    {
        std::cerr << "Input directory not specified\n";
        return usage(argv[0], 1);
    }
    if(sector.empty())
    {
        std::cerr << "Sector not specified\n";
        return usage(argv[0], 1);
    }
    if(outFileName.empty())
    {
        std::cerr << "Output file name not specified\n";
        return usage(argv[0], 1);
    }

    constexpr int wavelengths[] = {360,415,566,604,643,689};
    constexpr size_t numWLs = std::size(wavelengths);
    std::vector<std::vector<float>> dataPerWL(numWLs);
    ssize_t width = -1, height = -1, rowStride = -1;
    std::cerr << "Reading input data...\n";
    for(size_t wlN=0; wlN<std::size(wavelengths); ++wlN)
    {
        const auto filename = inDir + "/WAC_HAPKE_" + std::to_string(wavelengths[wlN])
                                    + "NM_" + sector + ".IMG";
        QFile file(filename.c_str());
        if(!file.open(QFile::ReadOnly))
        {
            std::cerr << "Failed to open " << filename << " for reading: "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
        QByteArray header(27360, '\0');
        if(file.read(header.data(), header.size()) != header.size())
        {
            std::cerr << "Failed to read header of " << filename << ": "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
        const auto currWidth = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
        const auto currHeight= QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
        if(currWidth==0 || currHeight==0)
        {
            std::cerr << "Failed to find image dimensions\n";
            return 1;
        }

        const auto currRowStride = currWidth;
        if(width<0)
            width = currWidth;
        if(height<0)
            height = currHeight;
        if(rowStride<0)
            rowStride = currRowStride;
        if(currWidth != width || currHeight != height || currRowStride != rowStride)
        {
            std::cerr << "Image dimensions or layout mismatch at " << filename << "\n";
            return 1;
        }

        dataPerWL[wlN].resize(rowStride*height);
        const qint64 sizeToRead = dataPerWL[wlN].size() * sizeof dataPerWL[wlN][0];
        if(file.read(reinterpret_cast<char*>(dataPerWL[wlN].data()), sizeToRead) != sizeToRead)
        {
            std::cerr << "Failed to read " << filename << ": "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
    }

    // These values were computed from CIE 1931 functions using triangular indicator functions
    // multiplied with solar irradiance at near-Earth orbit (result should be equivalent to using
    // linearly interpolated spectra of lunar albedo).
    constexpr double rgbsPerWL[][3] = {{0.123544388739659, -0.106861878682004, 0.752259001205208 },
                                       {-19.8366079265001,  27.8275045114262,  151.197577413019  },
                                       { 4.43497668468578,  157.373608612933,  37.079183527029   },
                                       { 145.775550325006,  11.7192895881426, -4.38388154872636  },
                                       { 73.8247972448887, -5.47353897213371, -0.811862808601448 },
                                       { 8.24022224911719, -0.85351964335178, -0.0605466729302762}};

    std::cerr << "Converting data to colors...\n";
    static_assert(std::size(rgbsPerWL) == numWLs);
    constexpr int bytesPerColor = 3;
    std::vector<double> out(rowStride*height*bytesPerColor);
    for(size_t n = 0; n < dataPerWL[0].size(); ++n)
    {
        double r = 0, g = 0, b = 0;
        bool good = true;
        for(size_t wlN = 0; wlN < numWLs; ++wlN)
        {
            const auto value = dataPerWL[wlN][n];
            if(value <= 0 || !std::isfinite(value))
            {
                good = false;
                break;
            }
            r += value * rgbsPerWL[wlN][0];
            g += value * rgbsPerWL[wlN][1];
            b += value * rgbsPerWL[wlN][2];
        }

        if(good)
        {
            double min = INFINITY, max = -INFINITY;
            for(size_t wlN = 0; wlN < numWLs; ++wlN)
            {
                const auto value = dataPerWL[wlN][n];
                if(value > max) max = value;
                if(value < min) min = value;
            }
            if(max / min > 5)
                good = false;
        }

        if(good)
        {
            if(dataPerWL[0][n] > dataPerWL[1][n])
                good=false;
            if(dataPerWL[numWLs-2][n] > dataPerWL[numWLs-1][n])
                good=false;
        }

        if(good)
        {
            for(size_t wlN = 1; wlN < numWLs - 1; ++wlN)
            {
                const double c = 1.00;
                if(dataPerWL[wlN-1][n] > c*dataPerWL[wlN][n] && c*dataPerWL[wlN][n] < dataPerWL[wlN+1][n])
                {
                    good=false;
                    break;
                }
                if(dataPerWL[wlN-1][n]*c < dataPerWL[wlN][n] && dataPerWL[wlN][n] > c*dataPerWL[wlN+1][n])
                {
                    good=false;
                    break;
                }
            }
        }

        if(good || markBadMode == MarkBadMode::None)
        {
            out[bytesPerColor*n+0] = r;
            out[bytesPerColor*n+1] = g;
            out[bytesPerColor*n+2] = b;
        }
        else
        {
            // Mark with magenta
            out[bytesPerColor*n+0] = 1/valueScale;
            out[bytesPerColor*n+1] = 0;
            out[bytesPerColor*n+2] = 1/valueScale;
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

    std::vector<uint8_t> outImgData(rowStride*height*3);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(std::abs(out[n]) * valueScale, 0., 1.));
    const QImage outImg(outImgData.data(), width, height, rowStride*3, QImage::Format_RGB888);
    std::cerr << "Saving output file with dimensions " << width << " × " << height << "... ";
    if(!outImg.save(outFileName.c_str()))
    {
        std::cerr << "Failed to save output file " << outFileName << "\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
