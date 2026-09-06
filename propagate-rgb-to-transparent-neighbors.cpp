#include <iostream>
#include <QImage>

int usage(const char*const argv0, const int ret)
{
    std::cerr << "Usage: " << argv0 << "[options...] inputImage outputImage\n"
                 "This tool locates transparent pixels neighboring opaque ones and propagates the\n"
                 "RGB values of the opaque pixels into the transparent ones. The goal is to avoid\n"
                 "a dark line appearing when interpolating a non-alpha-premultiplied image.\n"
                 "Options:\n"
                 "   -h, --help                         This help message\n"
                 "   -t, --threshold N                  Max alpha to consider transparent, [0..255]\n"
                 ;
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","0",false);

    QString inFileName;
    QString outFileName;
    int alphaThreshold = 0;

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
#define GO_TO_PARAM()                                                       \
            ++n;                                                            \
            if(n == argc)                                                   \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return 1;                                                   \
            }                                                               \
            do{}while(0)

        if(arg == "-h" || arg == "--help")
            return usage(argv[0], 0);
        else if(arg == "-t" || arg == "--threshold")
        {
            GO_TO_PARAM();
            alphaThreshold = std::stoul(argv[n]);
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(totalPositionalArgumentsFound < 2)
    {
        std::cerr << "Too few positional arguments supplied\n";
        return usage(argv[0], 1);
    }

    QImage in(inFileName);
    if(in.isNull())
    {
        std::cerr << "Failed to open input file\n";
        return 1;
    }
    in = in.convertToFormat(QImage::Format_RGBA8888);

    const auto inputData = in.bits();
    const auto inputWidth  = in.width();
    const auto inputHeight = in.height();
    const auto inputStrideInBytes = in.bytesPerLine();
    if(inputStrideInBytes % sizeof inputData[0])
    {
        std::cerr << "Row stride of " << inputStrideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto inputRowStride = inputStrideInBytes / sizeof inputData[0];

    auto out = in;
    const auto outputData = out.bits();

    ssize_t pixelsChanged = 0;
    ssize_t transparentPixelsFound = 0;
    for(ssize_t j = 0; j < inputHeight; ++j)
    {
        for(ssize_t i = 0; i < inputWidth; ++i)
        {
            const auto centerAlpha = inputData[j * inputRowStride + i * 4 + 3];
            if(centerAlpha > alphaThreshold) continue;
            ++transparentPixelsFound;
            double red = 0, green = 0, blue = 0, alpha = 0;
            int numPoints = 0;
            const int offsets[][2] = {{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
            for(const auto& [di,dj] : offsets)
            {
                const auto I = i + di, J = j + dj;
                if(I < 0 || I >= inputWidth || J < 0 || J >= inputHeight)
                    continue;
                red   += inputData[J * inputRowStride + I * 4 + 0] / 255.;
                green += inputData[J * inputRowStride + I * 4 + 1] / 255.;
                blue  += inputData[J * inputRowStride + I * 4 + 2] / 255.;
                alpha += inputData[J * inputRowStride + I * 4 + 3] / 255.;
                ++numPoints;
            }
            if(alpha == 0) continue;

            outputData[j * inputRowStride + i * 4 + 0] = std::clamp(red   / alpha, 0., 1.) * 255;
            outputData[j * inputRowStride + i * 4 + 1] = std::clamp(green / alpha, 0., 1.) * 255;
            outputData[j * inputRowStride + i * 4 + 2] = std::clamp(blue  / alpha, 0., 1.) * 255;

            ++pixelsChanged;
        }
    }

    std::cerr << "Altered " << pixelsChanged << " pixels out of " << transparentPixelsFound << " transparent ones\n";

    if(!out.save(outFileName))
    {
        std::cerr << "Failed to save output file " << outFileName.toStdString() << "\n";
    }
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}

