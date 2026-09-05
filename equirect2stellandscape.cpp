#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>

#include <thread>
#include <QImage>
#include <QVector2D>

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height,
             const size_t rowStride,
             const size_t subpixelIndex, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);

    const auto rawValue = data[x*4 + subpixelIndex + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height,
              const size_t rowStride, const size_t subpixelIndex,
              const double longitude, const double latitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-M_PI/2 <= latitude && latitude <= M_PI);

    const auto deltaLon = 2.*M_PI/width;
    const auto firstLon = (1.-width)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto roundedX = std::lround(x);

    const auto deltaLat = -M_PI/height;
    const auto firstLat = (1.-height)/2. * deltaLat;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto roundedY = std::lround(y);

    return fetch(data, width, height, rowStride, subpixelIndex, roundedX, roundedY);
}

template<typename T>
int findEmptyTopHeight(T const* data, const ssize_t width, const ssize_t height,
                       const size_t rowStride)
{
    ssize_t jMax = -1;
    for(ssize_t j = 0; j < height; ++j)
    {
        bool allEmpty = true;
        for(ssize_t i = 0; i < width; ++i)
        {
            const auto alpha = fetch(data, width, height, rowStride, 3, i, j);
            if(alpha != 0)
            {
                allEmpty = false;
                break;
            }
        }
        if(allEmpty && j > jMax)
            jMax = j;
    }
    const int totalEmptyLinesAtTop = jMax+1;
    return totalEmptyLinesAtTop;
}

size_t roundToNextPowerOfTwo(const size_t x)
{
    size_t p=1;
    while(p < x)
        p *= 2;
    return p;
}

int usage(const char*const argv0, const int ret)
{
    std::cerr << "Usage: " << argv0 << "[options...] inputImage outputDir\n"
                 "Options:\n"
                 "   -h, --help                         This help message\n"
                 "   -w, --width N                      Width in pixels of each side (default: 512)\n"
                 "   -n, --sides N                      Number of sides (default: 4)\n"
                 "   -s, --supersampling N              Number of supersampling samples (default: 1)\n"
                 "   -b, --side-bottom-angle-shift R    The value of decor_angle_shift (default: -30)\n"
                 "   --no-force-height-pot              Don't force side height to be a power-of-two\n"
                 "   --align-horizon                    Align horizon to the pixel grid (adjusts the\n"
                 "                                       value of -b)\n"
                 ;
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","0",false);

    QString inFileName;
    QString outDir;
    ssize_t sideWidth = 512;
    int numSides = 4;
    int numSamples = 1;
    double sideBottomAngularShift = -30;
    double startingLongitude = M_PI;
    bool forceHeightPOT = true;
    bool alignHorizonToPixelGrid = false;

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
                outDir = argv[n];
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
        else if(arg == "-w" || arg == "--width")
        {
            GO_TO_PARAM();
            sideWidth = std::stoul(argv[n]);
        }
        else if(arg == "-n" || arg == "--sides")
        {
            GO_TO_PARAM();
            numSides = std::stoul(argv[n]);
        }
        else if(arg == "-s" || arg == "--supersampling")
        {
            GO_TO_PARAM();
            numSamples = std::stoul(argv[n]);
        }
        else if(arg == "-b" || arg == "--side-bottom-angle-shift")
        {
            GO_TO_PARAM();
            sideBottomAngularShift = std::stod(argv[n]) * (M_PI / 180);
        }
        else if(arg == "--no-force-height-pot")
        {
            forceHeightPOT = false;
        }
        else if(arg == "--align-horizon")
        {
            alignHorizonToPixelGrid = true;
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

    const double pixelsPerRadianAtHorizon = sideWidth * numSides / (2*M_PI);
    // Let it overlap by about overlapGradientHeight pixels with the sides
    constexpr int overlapGradientHeight = 10;
    const double groundAngleShift = sideBottomAngularShift + overlapGradientHeight / pixelsPerRadianAtHorizon;

    const int emptyTopHeight = findEmptyTopHeight(inputData, inputWidth, inputHeight, inputRowStride);
    std::cerr << "Input map has " << emptyTopHeight << "px out of " << inputHeight
              << " empty (" << 100.*emptyTopHeight/inputHeight << "%) at the top\n";
    const double emptyTopAngularExtent = M_PI * emptyTopHeight / inputHeight;
    double sideAngularHeight = M_PI/2 - sideBottomAngularShift - emptyTopAngularExtent;
    int sideHeight = sideAngularHeight * pixelsPerRadianAtHorizon;
    if(forceHeightPOT)
        sideHeight = roundToNextPowerOfTwo(sideHeight);
    sideAngularHeight = sideHeight / pixelsPerRadianAtHorizon; // update using the rounded value of sideHeight

    if(alignHorizonToPixelGrid)
    {
        constexpr double maxSuperSampleShift = 0.5;
        constexpr double latitudeToAlign = 0;
        const double jMaxUnaligned = (sideBottomAngularShift - latitudeToAlign) / sideAngularHeight * sideHeight + sideHeight - 0.5 + maxSuperSampleShift;
        // For alignment the maximum value of j corresponding to latitudeToAlign must be a half-integer, because it's right between two pixels
        const double jMaxAligned = std::round(jMaxUnaligned + 0.5) - 0.5;
        sideBottomAngularShift = latitudeToAlign - (sideHeight - 0.5 - jMaxAligned) / sideHeight * sideAngularHeight;
    }

    std::vector<QVector2D> samples(numSamples);
    if(numSamples<=1)
        samples = {QVector2D(0,0)};
    else
    {
        std::mt19937 mt(std::random_device{}());
        std::uniform_real_distribution<double> dist(-0.5,0.5);
        static const auto rng=[&](){return dist(mt);};

        for(auto& sample : samples)
            sample = QVector2D(rng(), rng());
    }

    const double groundVShift = std::tan(groundAngleShift);
    const double sideAngularWidth = 2*M_PI / numSides;
    for(int side = -1; side < numSides; ++side)
    {
        const bool isGround = side == -1;
        std::cerr << "Processing " << (isGround ? "ground" : "side "+std::to_string(side)) << "...\n";
        const int currentSideHeight = isGround ? sideWidth : sideHeight;
        const auto dataSize = sideWidth*currentSideHeight*4;
        std::vector<double> outData(dataSize);
        for(ssize_t i = 0; i < sideWidth; ++i)
        {
            for(ssize_t j = 0; j < currentSideHeight; ++j)
            {
                const auto pixelPosInData = (i + j*sideWidth)*4;
                double red = 0, green = 0, blue = 0, alpha = 0;
                for(int sampleN = 0; sampleN < numSamples; ++sampleN)
                {
                    const double I = i + samples[sampleN].x();
                    const double x = 2 * ((I + 0.5) / sideWidth - 0.5); // range: (-1,1)
                    double longitude;
                    const double J = j + samples[sampleN].y();
                    const double y = 2 * ((J + 0.5) / currentSideHeight - 0.5); // range: (-1,1)
                    double latitude;
                    if(isGround)
                    {
                        const double z = groundVShift;
                        latitude = std::atan(z / std::sqrt(x*x+y*y));
                        longitude = std::atan2(y, x) - M_PI / 2;
                        if(longitude < -M_PI)
                            longitude += 2*M_PI;
                    }
                    else
                    {
                        latitude  = sideBottomAngularShift + (currentSideHeight - 0.5 - J) / currentSideHeight * sideAngularHeight;
                        longitude = startingLongitude + (side + (I + 0.5) / sideWidth) * sideAngularWidth;
                    }

                    red   += sample(inputData, inputWidth, inputHeight, inputRowStride, 0, longitude, latitude);
                    green += sample(inputData, inputWidth, inputHeight, inputRowStride, 1, longitude, latitude);
                    blue  += sample(inputData, inputWidth, inputHeight, inputRowStride, 2, longitude, latitude);
                    alpha += sample(inputData, inputWidth, inputHeight, inputRowStride, 3, longitude, latitude);
                }
                if(alpha > 0)
                {
                    outData[pixelPosInData + 0] = red   / alpha;
                    outData[pixelPosInData + 1] = green / alpha;
                    outData[pixelPosInData + 2] = blue  / alpha;
                    outData[pixelPosInData + 3] = alpha / numSamples;
                }
                else
                {
                    outData[pixelPosInData + 0] = red   / numSamples;
                    outData[pixelPosInData + 1] = green / numSamples;
                    outData[pixelPosInData + 2] = blue  / numSamples;
                    outData[pixelPosInData + 3] = alpha / numSamples;
                }
            }
        }
        if(!isGround)
        {
            // Fade out the bottom into the ground to avoid seams due to resolution or grid orientation mismatch
            for(ssize_t j = currentSideHeight - 1 - overlapGradientHeight; j < currentSideHeight; ++j)
            {
                const double alphaFactor = double(currentSideHeight - j) / (overlapGradientHeight + 1);
                for(ssize_t i = 0; i < sideWidth; ++i)
                {
                    const auto pixelPosInData = (i + j*sideWidth)*4;
                    outData[pixelPosInData + 3] *= alphaFactor;
                }
            }
        }


        using OutType = uchar;
        std::vector<OutType> outBits;
        outBits.reserve(outData.size());
        if(std::is_integral_v<OutType>)
        {
            for(auto v : outData)
                outBits.push_back(std::clamp(v,0.,1.)*std::numeric_limits<OutType>::max());
        }
        else
        {
            for(auto v : outData)
                outBits.push_back(v*std::numeric_limits<OutType>::max());
        }
        const QImage out(outBits.data(), sideWidth, currentSideHeight, sideWidth*4, in.format());
        const auto baseName = isGround ? "ground" : QString("/side%1").arg(side);
        const auto fileName = outDir + "/" + baseName + ".png";
        if(!out.save(fileName))
            std::cerr << "Failed to save output file " << fileName.toStdString() << "\n";
    }

    std::cerr << "Landscape parameters for the [landscape] section:\n\n";
    std::cout << "type = old_style\n";
    std::cout << "calibrated = true\n";
    std::cout << "nbsidetex = " << numSides << "\n";
    for(int side = 0; side < numSides; ++side)
        std::cout << "tex" << side << " = side" << side << ".png\n";
    std::cout << "nbside = " << numSides << "\n";
    for(int side = 0; side < numSides; ++side)
        std::cout << "side" << side << " = tex" << side << ":0:0:1:1\n";
    std::cout << "groundtex = ground.png\n";
    std::cout << "ground = groundtex:0:0:1:1\n";
    std::cout << "nb_decor_repeat = 1\n";
    std::cout << "decor_alt_angle = " << sideAngularHeight * (180 / M_PI) << "\n";
    std::cout << "decor_angle_shift = " << sideBottomAngularShift * (180 / M_PI) << "\n";
    std::cout << "decor_angle_rotatez = 0\n";
    std::cout << "ground_angle_shift = " << groundAngleShift * (180 / M_PI) << "\n";
    std::cout << "ground_angle_rotatez = 0\n";
    std::cout << "draw_ground_first = true\n";

}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}

