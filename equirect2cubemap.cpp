#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>

#include <thread>
#include <QImage>
#include <QVector2D>

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height,
             const size_t rowStride, const size_t channelsPerPixel,
             const size_t subpixelIndex, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);

    const auto rawValue = data[x*channelsPerPixel + subpixelIndex + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height,
              const size_t rowStride, const size_t channelsPerPixel,
              const size_t subpixelIndex, const double longitude, const double latitude)
{
    assert(-M_PI <= longitude && longitude <= M_PI);
    assert(-M_PI/2 <= latitude && latitude <= M_PI);

    const auto deltaLon = 2.*M_PI/width;
    const auto firstLon = (1.-width)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -M_PI/height;
    const auto firstLat = (1.-height)/2. * deltaLat;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

int usage(const char*const argv0, const int ret)
{
    std::cerr << "Usage: " << argv0 << "[options...] inputFile outputFilePattern cubeMapSideInPixels [supersamplingLevel]\n"
                 "Output file pattern must contain a %1 placeholder that will be replaced with face name when saving the faces\n"
                 "Options:\n"
                 "   -r                             Repeat each edge in both of its neighboring faces\n"
                 "   -e, --equiang, --equiangular   Make an equi-angular cube map\n"
                 ;
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    QString inFileName;
    QString outFileNamePattern;
    bool repeatEdgeInBothNeighbors = false;
    unsigned cubeMapSide = 0;
    unsigned numSamples = 1;
    bool equiangular = false;

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
                outFileNamePattern = argv[n];
                break;
            case 2:
                cubeMapSide = std::stoul(argv[n]);
                break;
            case 3:
                numSamples = std::stoul(argv[n]);
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
        else if(arg == "-r")
            repeatEdgeInBothNeighbors = true;
        else if(arg == "-e" || arg == "--equiang" || arg == "--equiangular")
            equiangular = true;
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(totalPositionalArgumentsFound < 3)
    {
        std::cerr << "Too few positional arguments supplied\n";
        return usage(argv[0], 1);
    }

    if(!outFileNamePattern.contains("%1"))
    {
        std::cerr << "File name pattern must contain a %1 placeholder that will be replaced with face name when saving the faces\n";
        return 1;
    }

    QImage in(inFileName);
    if(in.isNull())
    {
        std::cerr << "Failed to open input file\n";
        return 1;
    }
    in = in.convertToFormat(in.isGrayscale() ? QImage::Format_Grayscale8 : in.hasAlphaChannel() ? QImage::Format_ARGB32 : QImage::Format_RGB32);

    const auto data = in.bits();
    const auto width  = in.width();
    const auto height = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof data[0])
    {
        std::cerr << "Row stride of " << strideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto bitDepth = in.depth();
    assert(bitDepth % (8 * sizeof data[0]) == 0);
    const int channelsPerPixel = bitDepth / (8 * sizeof data[0]);
    const auto rowStride = strideInBytes / sizeof data[0];

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

#define FACE_LOOP(DATA_NAME, UV_ORIGIN, X_EXPR, Y_EXPR, Z_EXPR)                     \
    for(unsigned sampleN = 0; sampleN < numSamples; ++sampleN)                      \
    {                                                                               \
        for(size_t i = 0; i < cubeMapSide; ++i)                                     \
        {                                                                           \
            const auto I = i + samples[sampleN].x();                                \
            auto u = (UV_ORIGIN + I) / (cubeMapSide+2*UV_ORIGIN-1);                 \
            if(equiangular) u = (std::tan(M_PI/2*(u-0.5))+1)*0.5;                   \
            for(size_t j = 0; j < cubeMapSide; ++j)                                 \
            {                                                                       \
                const auto J = j + samples[sampleN].y();                            \
                auto v = (UV_ORIGIN+(cubeMapSide - 1 - J)) / (cubeMapSide+2*UV_ORIGIN-1); \
                if(equiangular) v = (std::tan(M_PI/2*(v-0.5))+1)*0.5;               \
                const auto pixelPosInData = (i + j*cubeMapSide)*channelsPerPixel;   \
                const auto x = X_EXPR;                                              \
                const auto y = Y_EXPR;                                              \
                const auto z = Z_EXPR;                                              \
                                                                                    \
                const auto longitude = atan2(y,x);                                  \
                const auto latitude = std::asin(z / std::sqrt(x*x+y*y+z*z));        \
                                                                                    \
                for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)   \
                {                                                                   \
                    DATA_NAME[pixelPosInData + subpixelN] +=                        \
                        sample(data, width, height, rowStride, channelsPerPixel,    \
                               subpixelN, longitude, latitude);                     \
                }                                                                   \
            }                                                                       \
        }                                                                           \
    }                                                                               \
    if(numSamples > 1)                                                              \
    {                                                                               \
        for(auto& v : DATA_NAME)                                                    \
            v *= 1./numSamples;                                                     \
    }                                                                               \
    do{}while(false)

    const auto dataSize = cubeMapSide*cubeMapSide*channelsPerPixel;
    std::vector<double> outDataNorth (dataSize);
    std::vector<double> outDataSouth (dataSize);
    std::vector<double> outDataWest  (dataSize);
    std::vector<double> outDataLon0  (dataSize);
    std::vector<double> outDataEast  (dataSize);
    std::vector<double> outDataLon180(dataSize);
    if(repeatEdgeInBothNeighbors)
    {
        std::thread n([&]{ FACE_LOOP(outDataNorth , 0,   u*2-1 ,   v*2-1 ,   1  ); });
        std::thread s([&]{ FACE_LOOP(outDataSouth , 0,   u*2-1 , -(v*2-1),  -1  ); });
        std::thread w([&]{ FACE_LOOP(outDataWest  , 0,   u*2-1 ,    -1   , v*2-1); });
        std::thread z([&]{ FACE_LOOP(outDataLon0  , 0,     1   ,   u*2-1 , v*2-1); });
        std::thread e([&]{ FACE_LOOP(outDataEast  , 0, -(u*2-1),     1   , v*2-1); });
        std::thread o([&]{ FACE_LOOP(outDataLon180, 0,    -1   , -(u*2-1), v*2-1); });
        n.join();
        s.join();
        w.join();
        z.join();
        e.join();
        o.join();
        // Our supersampling works in different areas on the edges, which yields
        // non-identical results in neighboring faces. This is not correct, but fixing
        // it seems a bit elaborate, while the improvement seems not too significant.
        // Instead, we'll copy the edge from one neighbor to another to make sure that
        // they are the same.
#define COPY_VERT_LINE_RIGHT_TO_LEFT(LEFT_DATA,RIGHT_DATA)                                  \
        for(unsigned n = 0; n < cubeMapSide; ++n)                                           \
        {                                                                                   \
            const auto pixelPosInSrc = channelsPerPixel*(n*cubeMapSide);                    \
            const auto pixelPosInDst = channelsPerPixel*(n*cubeMapSide+cubeMapSide-1);      \
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)               \
                LEFT_DATA[pixelPosInDst+subpixelN] = RIGHT_DATA[pixelPosInSrc+subpixelN];   \
        } do{}while(0)
        COPY_VERT_LINE_RIGHT_TO_LEFT(outDataLon0,outDataEast);
        COPY_VERT_LINE_RIGHT_TO_LEFT(outDataEast,outDataLon180);
        COPY_VERT_LINE_RIGHT_TO_LEFT(outDataLon180,outDataWest);
        COPY_VERT_LINE_RIGHT_TO_LEFT(outDataWest,outDataLon0);
        // North to sides: west
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)+n);
            const auto pixelPosInDst = channelsPerPixel*n;
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataWest[pixelPosInDst+subpixelN] = outDataNorth[pixelPosInSrc+subpixelN];
        }
        // North to sides: east
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide-1-n);
            const auto pixelPosInDst = channelsPerPixel*n;
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataEast[pixelPosInDst+subpixelN] = outDataNorth[pixelPosInSrc+subpixelN];
        }
        // North to sides: longitude 0
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide*cubeMapSide-1-n*cubeMapSide);
            const auto pixelPosInDst = channelsPerPixel*n;
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataLon0[pixelPosInDst+subpixelN] = outDataNorth[pixelPosInSrc+subpixelN];
        }
        // North to sides: longitude 180
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(n*cubeMapSide);
            const auto pixelPosInDst = channelsPerPixel*n;
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataLon180[pixelPosInDst+subpixelN] = outDataNorth[pixelPosInSrc+subpixelN];
        }
        // South to sides: west
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*n;
            const auto pixelPosInDst = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)+n);
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataWest[pixelPosInDst+subpixelN] = outDataSouth[pixelPosInSrc+subpixelN];
        }
        // South to sides: east
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide*cubeMapSide-1-n);
            const auto pixelPosInDst = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)+n);
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataEast[pixelPosInDst+subpixelN] = outDataSouth[pixelPosInSrc+subpixelN];
        }
        // South to sides: longitude 0
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide-1+n*cubeMapSide);
            const auto pixelPosInDst = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)+n);
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataLon0[pixelPosInDst+subpixelN] = outDataSouth[pixelPosInSrc+subpixelN];
        }
        // South to sides: longitude 180
        for(unsigned n = 0; n < cubeMapSide; ++n)
        {
            const auto pixelPosInSrc = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)-n*cubeMapSide);
            const auto pixelPosInDst = channelsPerPixel*(cubeMapSide*(cubeMapSide-1)+n);
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
                outDataLon180[pixelPosInDst+subpixelN] = outDataSouth[pixelPosInSrc+subpixelN];
        }
    }
    else
    {
        std::thread n([&]{ FACE_LOOP(outDataNorth , 0.5,   u*2-1 ,   v*2-1 ,   1  ); });
        std::thread s([&]{ FACE_LOOP(outDataSouth , 0.5,   u*2-1 , -(v*2-1),  -1  ); });
        std::thread w([&]{ FACE_LOOP(outDataWest  , 0.5,   u*2-1 ,    -1   , v*2-1); });
        std::thread z([&]{ FACE_LOOP(outDataLon0  , 0.5,     1   ,   u*2-1 , v*2-1); });
        std::thread e([&]{ FACE_LOOP(outDataEast  , 0.5, -(u*2-1),     1   , v*2-1); });
        std::thread o([&]{ FACE_LOOP(outDataLon180, 0.5,    -1   , -(u*2-1), v*2-1); });
        n.join();
        s.join();
        w.join();
        z.join();
        e.join();
        o.join();
    }

    const QString faceTypes[]={"3-west", "0-lon0", "2-east", "1-lon180", "5-north","4-south"};
    unsigned faceN=0;
    for(const auto& outData : {outDataWest, outDataLon0, outDataEast, outDataLon180, outDataNorth, outDataSouth})
    {
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
        const QImage out(outBits.data(), cubeMapSide, cubeMapSide, cubeMapSide*channelsPerPixel, in.format());
        const auto fileName = outFileNamePattern.arg(faceTypes[faceN]);
        if(!out.save(fileName))
        {
            std::cerr << "Failed to save output file " << fileName.toStdString() << "\n";
            return 1;
        }
        ++faceN;
    }
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
