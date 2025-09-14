#include <iostream>
#include <stdexcept>
#include <QDebug>
#include <QImage>
#include <glm/glm.hpp>

uint64_t toUInt64(glm::vec3 value)
{
    const uint64_t r = std::clamp(int(value[0]), 0, 65535);
    const uint64_t g = std::clamp(int(value[1]), 0, 65535);
    const uint64_t b = std::clamp(int(value[2]), 0, 65535);
    return 0xffffULL << 48 | b << 32 | g << 16 | r;
}

glm::vec3 fetch(const uint64_t*const data, const ssize_t stride,
                const ssize_t width, const ssize_t height,
                ssize_t i, ssize_t j)
{
    // Wrap around for longitude
    while(i < 0) i += width;
    while(i >= width) i -= width;
    // Clamping for latitude
    if(j < 0) j = 0;
    if(j >= height) j = height - 1;

    const auto rgbx = data[j*stride + i];
    return {float(rgbx&0xffff), float((rgbx >> 16) & 0xffff), float((rgbx >> 32) & 0xffff)};
}

glm::vec3 sample(const uint64_t*const data, const ssize_t stride,
                 const ssize_t width, const ssize_t height,
                 const double i, const double j)
{
    const auto iLeft = std::floor(i), jTop = std::floor(j);
    const float alphaI = i - iLeft, alphaJ = j - jTop;
    const auto vTopLeft     = fetch(data, stride, width, height, iLeft, jTop);
    const auto vBottomLeft  = fetch(data, stride, width, height, iLeft, jTop + 1);
    const auto vTopRight    = fetch(data, stride, width, height, iLeft + 1, jTop);
    const auto vBottomRight = fetch(data, stride, width, height, iLeft + 1, jTop + 1);
    const auto top = vTopLeft*(1 - alphaI) + vTopRight*alphaI;
    const auto bottom = vBottomLeft*(1 - alphaI) + vBottomRight*alphaI;
    return top*(1 - alphaJ) + bottom*alphaJ;
}

std::vector<glm::ivec2> parseSeparatingLines(const QString& spec)
{
    const auto coords = spec.split(",");
    if(coords.size() % 2)
        throw std::runtime_error("Number of coordinates in the polyline must be even: N entries of (x,y)");
    std::vector<glm::ivec2> out;
    for(int n = 0; n < coords.size(); n += 2)
    {
        bool ok1 = false, ok2 = false;
        out.emplace_back(coords[n].toInt(&ok1), coords[n+1].toInt(&ok2));
        if(!ok1 || !ok2)
            throw std::runtime_error(("Failed to parse line coordinates "+
                                      coords[n]+" and "+coords[n+1]).toStdString());
        if(n > 0 && out[out.size() - 2][0] >= out.back()[0])
            throw std::runtime_error("Line x coordinates aren't monotonic at point "+std::to_string(n/2));
        if(out.back()[0] < 0)
            throw std::runtime_error("Negative x coordinate of a line point detected");
    }
    return out;
}

int usage(const char*const argv0, const int ret)
{
        std::cerr << "Usage: " << argv0 << 1+R"(
 {options...} map-left-or-top.png map-right-or-bottom.png combined-map.png
Options:
 --help                         Show this help message and quit
 --row NUM                      Connect the maps at row NUM
 --col NUM                      Connect the maps at column NUM
 --horiz-lines x1,y1,...,xN,yN  Connect the maps top to bottom at the polyline
Only one of --row, --col, --horiz-lines option can be used.
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","8192",false);

    if(argc == 1)
        return usage(argv[0],1);

    QString imageLeftFileName;
    QString imageRightFileName;
    QString outFileName;
    ssize_t connectingColumn = 0;
    ssize_t connectingRow = 0;
    std::vector<glm::ivec2> lines;
    int totalPositionalArgumentsFound = 0;

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
                imageLeftFileName = argv[n];
                break;
            case 1:
                imageRightFileName = argv[n];
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
        if(arg == "--help")
        {
            return usage(argv[0], 0);
        }
        else if(arg == "--col")
        {
            REQUIRE_PARAM();
            connectingColumn = std::stoul(argv[++n]);
        }
        else if(arg == "--row")
        {
            REQUIRE_PARAM();
            connectingRow = std::stoul(argv[++n]);
        }
        else if(arg == "--horiz-lines")
        {
            REQUIRE_PARAM();
            lines = parseSeparatingLines(argv[++n]);
        }
        else
        {
            std::cerr << "Unknown option " << arg << "\n";
            return usage(argv[0], 1);
        }
    }
#undef REQUIRE_PARAM
    if(imageLeftFileName.isEmpty() || imageRightFileName.isEmpty() || outFileName.isEmpty())
        return usage(argv[0], 1);
    if((connectingColumn > 0) + (connectingRow > 0) + !lines.empty() > 1)
    {
        std::cerr << "Only one of --row, --col, --horiz-lines option can be used\n";
        return 1;
    }
    if(connectingRow == 0 && connectingColumn == 0 && lines.empty())
    {
        std::cerr << "Neither --col nor --row nor --horiz-lines were specified\n";
        return 1;
    }

    std::cerr << "Reading first image... ";
    const auto imageLeft = QImage(imageLeftFileName).convertToFormat(QImage::Format_RGBX64);
    if(imageLeft.isNull())
    {
        std::cerr << "FAILED\n";
        return 1;
    }
    else
        std::cerr << "done\n";

    std::cerr << "Reading second image... ";
    const auto imageRight = QImage(imageRightFileName).convertToFormat(QImage::Format_RGBX64);
    if(imageRight.isNull())
    {
        std::cerr << "FAILED\n";
        return 1;
    }
    else
        std::cerr << "done\n";

    if(imageLeft.size() != imageRight.size())
    {
        std::cerr << "Error: input images have different sizes\n";
        return 1;
    }

    const ssize_t W = imageLeft.width(), H = imageLeft.height();
    QImage imageOut(W, H, QImage::Format_RGBX64);

    std::cerr << "Combining...\n";
    const auto dataLeft  = reinterpret_cast<const uint64_t*>(imageLeft.bits());
    const auto dataRight = reinterpret_cast<const uint64_t*>(imageRight.bits());
    const auto dataOut = reinterpret_cast<uint64_t*>(imageOut.bits());
    const ssize_t strideLeft  = imageLeft.bytesPerLine()  / sizeof dataLeft[0];
    const ssize_t strideRight = imageRight.bytesPerLine() / sizeof dataRight[0];
    const ssize_t strideOut   = imageOut.bytesPerLine()    / sizeof dataOut[0];
    constexpr int ssFactor = 8;
    const ssize_t stripHalfWidth = 4;
    using std::sqrt;
    using std::abs;
    if(connectingColumn)
    {
        ssize_t rowIndexLeft = 0, rowIndexRight = 0, rowIndexOut = 0;
        for(ssize_t j = 0; j < H; ++j)
        {
            for(ssize_t i = 0; i < W; ++i)
            {
                if(i < connectingColumn - stripHalfWidth)
                {
                    dataOut[rowIndexOut + i] = dataLeft [rowIndexLeft + i];
                    continue;
                }
                else if(i > connectingColumn + stripHalfWidth)
                {
                    dataOut[rowIndexOut + i] = dataRight[rowIndexRight + i];
                    continue;
                }
                else
                {
                    const auto nearestLeft  = fetch(dataLeft,  strideLeft,  W, H, i, j);
                    const auto nearestRight = fetch(dataRight, strideRight, W, H, i, j);
                    glm::vec3 value = {0,0,0};
                    for(int ti = 0; ti < ssFactor; ++ti)
                    {
                        const float t = float(ti + 0.5f) / ssFactor;
                        for(int si = 0; si < ssFactor; ++si)
                        {
                            const float s = float(si + 0.5f) / ssFactor;
                            const auto alpha = std::clamp(1 - sqrt(abs(i + s - connectingColumn) / stripHalfWidth),
                                                          0.f, 1.f);
                            const auto sampLeft  = sample(dataLeft,  strideLeft,  W, H, i + s, j + t);
                            const auto sampRight = sample(dataRight, strideRight, W, H, i + s, j + t);
                            const auto samp = sampLeft[0] < sampRight[0] ? sampRight : sampLeft;
                            const auto partialVal = i + s < connectingColumn ? samp * alpha + nearestLeft  * (1 - alpha)
                                                                             : samp * alpha + nearestRight * (1 - alpha);
                            value += partialVal;
                        }
                    }
                    value = round(value / float(ssFactor*ssFactor));
                    dataOut[rowIndexOut + i] = toUInt64(value);
                }
            }
            rowIndexLeft += strideLeft;
            rowIndexRight += strideRight;
            rowIndexOut += strideOut;
        }
    }
    else if(connectingRow)
    {
        ssize_t rowIndexTop = 0, rowIndexBottom = 0, rowIndexOut = 0;
        auto& dataTop = dataLeft;
        auto& dataBottom = dataRight;
        auto& strideTop = strideLeft;
        auto& strideBottom = strideRight;
        for(ssize_t j = 0; j < H; ++j)
        {
            if(j < connectingRow - stripHalfWidth)
            {
                std::copy_n(dataTop + rowIndexTop, W, dataOut + rowIndexOut);
            }
            else if(j > connectingRow + stripHalfWidth)
            {
                std::copy_n(dataBottom + rowIndexBottom, W, dataOut + rowIndexOut);
            }
            else
            {
                for(ssize_t i = 0; i < W; ++i)
                {
                    const auto nearestTop    = fetch(dataTop,    strideTop,    W, H, i, j);
                    const auto nearestBottom = fetch(dataBottom, strideBottom, W, H, i, j);
                    glm::vec3 value = {0,0,0};
                    for(int ti = 0; ti < ssFactor; ++ti)
                    {
                        const float t = float(ti + 0.5f) / ssFactor;
                        for(int si = 0; si < ssFactor; ++si)
                        {
                            const float s = float(si + 0.5f) / ssFactor;
                            const auto alpha = std::clamp(1 - sqrt(abs(j + t - connectingRow) / stripHalfWidth),
                                                          0.f, 1.f);
                            const auto sampTop    = sample(dataTop,    strideTop,    W, H, i + s, j + t);
                            const auto sampBottom = sample(dataBottom, strideBottom, W, H, i + s, j + t);
                            const auto samp = sampTop[0] < sampBottom[0] ? sampBottom : sampTop;
                            const auto partialVal = j + t < connectingRow ? samp * alpha + nearestTop    * (1 - alpha)
                                                                          : samp * alpha + nearestBottom * (1 - alpha);
                            value += partialVal;
                        }
                    }
                    value = round(value / float(ssFactor*ssFactor));
                    dataOut[rowIndexOut + i] = toUInt64(value);
                }
            }
            rowIndexTop += strideTop;
            rowIndexBottom += strideBottom;
            rowIndexOut += strideOut;
        }
    }
    else
    {
        std::vector<float> ys;
        ys.reserve(W);
        ys.push_back(lines[0][1]);
        for(unsigned n = 1; n < lines.size(); ++n)
        {
            const auto p0 = glm::vec2(lines[n-1][0], lines[n-1][1]);
            const auto p1 = glm::vec2(lines[n][0], lines[n][1]);
            const auto L = p1[0] - p0[0];
            for(int x = lines[n-1][0]; x < lines[n][0]; ++x)
            {
                const auto t = (x - p0[0]) / L;
                const auto y = p0[1] + t * (p1[1] - p0[1]);
                ys.push_back(y);
            }
        }
        if(ssize_t(ys.size()) != W)
            throw std::logic_error("ys.size = "+std::to_string(ys.size())+" != W = "+std::to_string(W));

        const auto [minY, maxY] = std::minmax_element(ys.begin(), ys.end());

        ssize_t rowIndexTop = 0, rowIndexBottom = 0, rowIndexOut = 0;
        auto& dataTop = dataLeft;
        auto& dataBottom = dataRight;
        auto& strideTop = strideLeft;
        auto& strideBottom = strideRight;
        for(ssize_t j = 0; j < H; ++j)
        {
            if(j < *minY)
            {
                std::copy_n(dataTop + rowIndexTop, W, dataOut + rowIndexOut);
            }
            else if(j > *maxY)
            {
                std::copy_n(dataBottom + rowIndexBottom, W, dataOut + rowIndexOut);
            }
            else
            {
                for(ssize_t i = 0; i < W; ++i)
                {
                    // FIXME: no filtering between the images for now
                    if(j < ys[i])
                        dataOut[rowIndexOut + i] = dataTop[rowIndexTop + i];
                    else
                        dataOut[rowIndexOut + i] = dataBottom[rowIndexBottom + i];
                }
            }
            rowIndexTop += strideTop;
            rowIndexBottom += strideBottom;
            rowIndexOut += strideOut;
        }
    }

    std::cerr << "Saving output file... ";
    if(!imageOut.save(outFileName))
    {
        std::cerr << "Failed to save output file\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::logic_error const& ex)
{
    std::cerr << "Internal error: " << ex.what() << "\n";
    return 1;
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
