#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>

#include <QImage>
#include <QVector2D>

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height,
             const size_t rowStride, const size_t bytesPerPixel,
             const size_t subpixelIndex, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, 0, height-1);

    const auto rawValue = data[x*bytesPerPixel + subpixelIndex + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height,
              const size_t rowStride, const size_t bytesPerPixel,
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

    const auto pTopLeft     = fetch(data, width, height, rowStride, bytesPerPixel, subpixelIndex, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, bytesPerPixel, subpixelIndex, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, bytesPerPixel, subpixelIndex, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, bytesPerPixel, subpixelIndex, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

int main(int argc, char** argv)
try
{
    if(argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " inputFile outputFilePattern cubeMapSideInPixels [supersamplingLevel]\n";
        return 1;
    }

    const char*const inFileName  = argv[1];
    const auto outFileNamePattern = QString(argv[2]);
    if(!outFileNamePattern.contains("%1"))
    {
        std::cerr << "File name pattern must contain a %1 placeholder that will be replaced with face name when saving the faces\n";
        return 1;
    }

    const auto cubeMapSide = std::stoul(argv[3]);
    const unsigned numSamples = argc==4 ? 1 : std::stoul(argv[4]);

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
    if(bitDepth % 8)
    {
        std::cerr << bitDepth << " bit per pixel is not a multiple of octet, this is not supported\n";
        return 1;
    }
    const auto bytesPerPixel = bitDepth / 8;
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

#define FACE_LOOP(DATA_NAME, X_EXPR, Y_EXPR, Z_EXPR)                            \
    std::vector<double> DATA_NAME(cubeMapSide*cubeMapSide*bytesPerPixel);       \
    for(unsigned sampleN = 0; sampleN < numSamples; ++sampleN)                  \
    {                                                                           \
        for(size_t i = 0; i < cubeMapSide; ++i)                                 \
        {                                                                       \
            const auto I = i + samples[sampleN].x();                            \
            const auto u = (0.5 + I) / cubeMapSide;                             \
            for(size_t j = 0; j < cubeMapSide; ++j)                             \
            {                                                                   \
                const auto J = j + samples[sampleN].y();                        \
                const auto v = (0.5+(cubeMapSide - 1 - J)) / cubeMapSide;       \
                const auto pixelPosInData = (i + j*cubeMapSide)*bytesPerPixel;  \
                const auto x = X_EXPR;                                          \
                const auto y = Y_EXPR;                                          \
                const auto z = Z_EXPR;                                          \
                                                                                \
                const auto longitude = atan2(y,x);                              \
                const auto latitude = std::asin(z / std::sqrt(x*x+y*y+z*z));    \
                                                                                \
                for(int subpixelN = 0; subpixelN < bytesPerPixel; ++subpixelN)  \
                {                                                               \
                    DATA_NAME[pixelPosInData + subpixelN] +=                    \
                        sample(data, width, height, rowStride, bytesPerPixel,   \
                               subpixelN, longitude, latitude);                 \
                }                                                               \
            }                                                                   \
        }                                                                       \
    }                                                                           \
    if(numSamples > 1)                                                          \
    {                                                                           \
        for(auto& v : DATA_NAME)                                                \
            v *= 1./numSamples;                                                 \
    }                                                                           \
    do{}while(false)

    FACE_LOOP(outDataNorth ,   u*2-1 ,   v*2-1 ,   1);
    FACE_LOOP(outDataSouth ,   u*2-1 , -(v*2-1),  -1);
    FACE_LOOP(outDataWest  ,   u*2-1 ,    -1   , v*2-1);
    FACE_LOOP(outDataLon0  ,     1   ,   u*2-1 , v*2-1);
    FACE_LOOP(outDataEast  , -(u*2-1),     1   , v*2-1);
    FACE_LOOP(outDataLon180,    -1   ,   u*2-1 , v*2-1);

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
        const QImage out(outBits.data(), cubeMapSide, cubeMapSide, cubeMapSide*bytesPerPixel, in.format());
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
