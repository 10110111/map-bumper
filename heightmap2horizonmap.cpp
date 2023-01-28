#include <cmath>
#include <random>
#include <thread>
#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <QImage>
#include "timing.hpp"

template<typename T> auto sqr(T const& x) { return x*x; }

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height, const size_t rowStride,
             const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);

    return data[x + y*rowStride];
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

    const auto deltaLat = -M_PI/(height-1);
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
    if(argc != 6 && argc != 7)
    {
        std::cerr << "Usage: " << argv[0] << " sphereRadiusInKM kmPerUnitValue inputFile outputFileName outputHeightInPixels\n";
        return 1;
    }

    const auto sphereRadius = std::stod(argv[1]);
    const auto kmPerUnit = std::stod(argv[2]);
    const char*const inFileName  = argv[3];
    const char*const outFileName = argv[4];

    const auto outputHeight = std::stoul(argv[5]);
    const auto outputWidth = 2*outputHeight;

    QImage in(inFileName);
    if(in.isNull())
    {
        std::cerr << "Failed to open input file\n";
        return 1;
    }
    if(!in.isGrayscale())
    {
        std::cerr << "Input image is not grayscale. Can't interpret it as a height map.\n";
        return 1;
    }
    in = in.convertToFormat(QImage::Format_Grayscale16);

    const auto data = reinterpret_cast<uint16_t*>(in.bits());
    const auto inputWidth  = in.width();
    const auto inputHeight = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof data[0])
    {
        std::cerr << "Row stride of " << strideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto rowStride = strideInBytes / sizeof data[0];

    const auto maxDataValue = *std::max_element(data, data+inputWidth*inputHeight);
    const auto maxAltitude = kmPerUnit * maxDataValue;
    const auto maxRadiusSquared = sqr(sphereRadius+maxAltitude);

    const size_t outBytesPerPixel = 4;
    std::vector<double> outData(outputWidth*outputHeight*outBytesPerPixel);
    std::atomic<unsigned> rowsDone{0};
    auto work = [inputWidth,inputHeight,outputWidth,outputHeight,
                 sphereRadius,kmPerUnit,rowStride,maxRadiusSquared,
                 data,outData=outData.data(),&rowsDone](const size_t jMin, const size_t jMax)
    {
        auto time0 = std::chrono::steady_clock::now();
        size_t rowsDoneInThisThreadAfterLastUpdate = 0;
        for(size_t j = jMin; j < jMax; ++j, ++rowsDoneInThisThreadAfterLastUpdate)
        {
            const auto v = (0.5+(outputHeight - 1 - j)) / outputHeight;
            for(size_t i = 0; i < outputWidth; ++i)
            {
                const auto u = (0.5 + i) / outputWidth;
                const auto pixelPosInData = (i + j*outputWidth)*outBytesPerPixel;

                using namespace glm;
                using namespace std;

                const double raySourceLon = (2*u-1)*M_PI;
                const double raySourceLat  = (2*v-1)*(M_PI/2);

                const double deltaLon = (2*M_PI/inputWidth) / cos(raySourceLat);
                const double eastLon = raySourceLon + deltaLon;
                const double eastLat = raySourceLat;

                const auto deltaLat = M_PI/inputHeight;
                const double northLat = raySourceLat + deltaLat;
                const double northLon = raySourceLon;

                const double raySourceHeight = kmPerUnit*sample(data, inputWidth, inputHeight, rowStride, raySourceLon, raySourceLat);
                const double raySourceRadius = sphereRadius + raySourceHeight;

                const double eastHeight = kmPerUnit*sample(data, inputWidth, inputHeight, rowStride, eastLon, eastLat);
                const double eastRadius = sphereRadius + eastHeight;

                const double northHeight = kmPerUnit*sample(data, inputWidth, inputHeight, rowStride, northLon, northLat);
                const double northRadius = sphereRadius + northHeight;

                const dvec3 zenithAtRaySource= dvec3(cos(raySourceLon) * cos(raySourceLat),
                                                     sin(raySourceLon) * cos(raySourceLat),
                                                     sin(raySourceLat));
                const dvec3 raySourcePoint = raySourceRadius * zenithAtRaySource;

                const dvec3 eastPoint = eastRadius * dvec3(cos(eastLon) * cos(eastLat),
                                                           sin(eastLon) * cos(eastLat),
                                                           sin(eastLat));
                const dvec3 northPoint = northRadius * dvec3(cos(northLon) * cos(northLat),
                                                             sin(northLon) * cos(northLat),
                                                             sin(northLat));
                const dvec3 deltaEast = eastPoint - raySourcePoint;
                const dvec3 deltaNorth = northPoint - raySourcePoint;

                unsigned rayIndex = 0;
                for(const dvec3 rayDir : {deltaNorth, deltaEast, -deltaNorth, -deltaEast})
                {
                    double sinHorizonElevation = -M_PI/2;
                    for(dvec3 rayPoint = raySourcePoint + rayDir; dot(rayPoint, rayPoint) <= maxRadiusSquared; rayPoint += rayDir)
                    {
                        const dvec3 zenithAtRayPoint = normalize(rayPoint);

                        const double longitude = atan2(rayPoint.y, rayPoint.x);
                        const double latitude = asin(zenithAtRayPoint.z);
                        const double altitude = kmPerUnit*sample(data, inputWidth, inputHeight, rowStride, longitude, latitude);

                        const dvec3 pointAtAlt = zenithAtRayPoint*(sphereRadius+altitude);
                        const double sinElev = dot(normalize(pointAtAlt-raySourcePoint), zenithAtRaySource); 
                        if(sinElev > sinHorizonElevation)
                            sinHorizonElevation = sinElev;
                    }
                    outData[pixelPosInData + rayIndex] = (sign(sinHorizonElevation)*sqrt(abs(sinHorizonElevation)))/2+0.5;
                    ++rayIndex;
                }
            }
            auto time1 = std::chrono::steady_clock::now();
            if(time1 - time0 > std::chrono::seconds(5))
            {
                rowsDone += rowsDoneInThisThreadAfterLastUpdate;
                rowsDoneInThisThreadAfterLastUpdate = 0;
                std::cerr << std::round(rowsDone * 10000. / outputHeight)/100. << "% done\n";
                time0 = time1;
            }
        }
    };

    auto time0 = std::chrono::steady_clock::now();
    const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for(size_t n = 0; n < numThreads; ++n)
    {
        const size_t jMin = outputHeight / numThreads * n;
        const size_t jMax = n+1 < numThreads ? outputHeight / numThreads * (n+1) : outputHeight;
        threads.emplace_back(work, jMin, jMax);
    }
    for(auto& thread : threads)
        thread.join();
    auto time1 = std::chrono::steady_clock::now();
    std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";

    using OutType = uchar;
    std::vector<OutType> outBits;
    outBits.reserve(outData.size());
    for(auto v : outData)
        outBits.push_back(std::clamp(v,0.,1.)*std::numeric_limits<OutType>::max());

    const QImage out(outBits.data(), outputWidth, outputHeight, outputWidth*outBytesPerPixel, QImage::Format_RGBA8888);
    if(!out.save(outFileName))
    {
        std::cerr << "Failed to save output file " << outFileName << "\n";
        return 1;
    }
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}

