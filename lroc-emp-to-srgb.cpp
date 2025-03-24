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

const float* getSurroundingColor(const float*const data,
                                  const ssize_t width, const ssize_t height, const int bytesPerColor,
                                  const ssize_t i0, const ssize_t j0)
{
    struct Color
    {
        float r,g,b;
        Color(const float* c) : r(c[0]), g(c[1]), b(c[2]) {}
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
 --bad-to-black             Fill bad pixels with black
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
        Black,
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
        else if(arg == "--bad-to-black")
        {
            markBadMode = MarkBadMode::Black;
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
    if(sector.size() != 9)
    {
        std::cerr << "Sector should be 9 characters long, e.g. E350S0450\n";
        return 1;
    }
    if(outFileName.empty())
    {
        std::cerr << "Output file name not specified\n";
        return usage(argv[0], 1);
    }

    std::cerr << "Reading input data...\n";
    constexpr int wavelength = 643;
    const auto filename = inDir + "/WAC_EMP_" + std::to_string(wavelength)
                                + "NM_" + sector + "_304P.IMG";
    QFile file(filename.c_str());
    if(!file.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << filename << " for reading: "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }
    QByteArray header(1024, '\0');
    if(file.read(header.data(), header.size()) != header.size())
    {
        std::cerr << "Failed to read header of " << filename << ": "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }
    const auto headerSize = QRegularExpression("RECORD_BYTES\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    if(headerSize == 0)
    {
        std::cerr << "Failed to determine header size\n";
        return 1;
    }
    header = QByteArray(headerSize, '\0');
    file.seek(0);
    if(file.read(header.data(), header.size()) != header.size())
    {
        std::cerr << "Failed to read header of " << filename << ": "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }

    const ssize_t width = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    const ssize_t height= QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    if(width==0 || height==0)
    {
        std::cerr << "Failed to find image dimensions\n";
        return 1;
    }

    const auto rowStride = width;
    std::vector<float> data(rowStride*height);
    const qint64 sizeToRead = data.size() * sizeof data[0];
    if(file.read(reinterpret_cast<char*>(data.data()), sizeToRead) != sizeToRead)
    {
        std::cerr << "Failed to read " << filename << ": "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }

    std::cerr << "Converting data to colors...\n";
    std::vector<float> out(rowStride*height);
    for(size_t n = 0; n < data.size(); ++n)
    {
        bool good = true;
        const auto value = data[n];
        if(value <= 0 || !std::isfinite(value))
            good = false;
        if(good || markBadMode == MarkBadMode::None)
            out[n] = value;
        else if(markBadMode == MarkBadMode::Black)
            out[n] = 0;
        else
            out[n] = -1;
    }

    if(markBadMode == MarkBadMode::Surroundings)
    {
        std::cerr << "Replacing bad pixels with some surrounding colors...\n";
        for(ssize_t j = 0; j < height; ++j)
        {
            const ssize_t linePos = j*width;
            for(ssize_t i = 0; i < width; ++i)
            {
                const auto pos = (linePos + i);
                if(out[pos] <= 0)
                {
                    const auto*const surroundingColor = getSurroundingColor(out.data(), width, height, 1, i, j);
                    assert(surroundingColor);
                    if(!surroundingColor) continue;
                    // Mark fixed colors with a negative sign to prevent using them as surroundings.
                    // The negativity will be removed by abs() when saving the image.
                    out[pos] = -surroundingColor[0];
                }
            }
        }
        // Blur the replaced pixels a bit to better fit into the surroundings
        for(ssize_t j = 0; j < height; ++j)
        {
            const ssize_t linePos = j*width;
            for(ssize_t i = 0; i < width; ++i)
            {
                const auto pos = (linePos + i);
                if(out[pos] < 0)
                {
                    using std::abs;
                    const auto left = std::max(i-1, ssize_t(0));
                    const auto right = std::min(i+1, width-1);
                    const auto top = std::max(j-1, ssize_t(0));
                    const auto bottom = std::min(j+1, height-1);
                    const auto B = 1;
                    out[pos] = abs(out[(top*width+left)*B]) + abs(out[(top*width+i)*B]) + abs(out[(top*width+right)*B])+
                                 abs(out[((j+0)*width+left)*B]) + abs(out[((j+0)*width+i)*B]) + abs(out[((j+0)*width+right)*B])+
                                 abs(out[(bottom*width+left)*B]) + abs(out[(bottom*width+i)*B]) + abs(out[(bottom*width+right)*B]);
                    out[pos] /= 9;
                }
            }
        }
    }

    std::vector<uint8_t> outImgData(rowStride*height);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(std::abs(out[n]) * valueScale, 0., 1.));
    const QImage outImg(outImgData.data(), width, height, rowStride, QImage::Format_Grayscale8);
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
