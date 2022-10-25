#include <cmath>
#include <iostream>
#include <algorithm>

#include <QImage>

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height, const size_t rowStride,
             const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, 0, height-1);

    const auto rawValue = data[x + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height, const size_t rowStride,
              const double longitude, const double latitude)
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

    const auto pTopLeft     = fetch(data, width, height, rowStride, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

int main(int argc, char** argv)
try
{
    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " inputFile outputFile cubeMapSideInPixels\n";
        return 1;
    }

    const char*const inFileName  = argv[1];
    const char*const outFileName = argv[2];
    const auto cubeMapSide = std::stoul(argv[3]);

    QImage in(inFileName);
    if(in.isNull())
    {
        std::cerr << "Failed to open input file\n";
        return 1;
    }
    in = in.convertToFormat(QImage::Format_Grayscale8);

    const auto data = in.bits();
    const auto width  = in.width();
    const auto height = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof data[0])
    {
        std::cerr << "Row stride of " << strideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto rowStride = strideInBytes / sizeof data[0];

    std::vector<double> outData(cubeMapSide*cubeMapSide);
    for(size_t i = 0; i < cubeMapSide; ++i)
    {
        const auto u = (0.5+i)/cubeMapSide;
        for(size_t j = 0; j < cubeMapSide; ++j)
        {
            const auto v = (0.5+(cubeMapSide-1-j))/cubeMapSide;

            const auto x = u*2-1;
            const auto y = v*2-1;
            const auto z = 1;

            const auto longitude = atan2(y,x);
            const auto latitude = std::asin(z / std::sqrt(x*x+y*y+z*z));

            outData[i + j*cubeMapSide] = sample(data, width, height, rowStride, longitude, latitude);
        }
    }

    using OutType = uchar;
    std::vector<OutType> outBits;
    outBits.reserve(outData.size());
    for(auto v : outData)
        outBits.push_back(v*std::numeric_limits<OutType>::max());
    QImage out(outBits.data(), cubeMapSide, cubeMapSide, cubeMapSide, QImage::Format_Grayscale8);
    if(!out.save(outFileName))
    {
        std::cerr << "Failed to save output file\n";
        return 1;
    }
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
