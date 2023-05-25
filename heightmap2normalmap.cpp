#include <cmath>
#include <random>
#include <atomic>
#include <thread>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <QImage>
#include "timing.hpp"

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

//! Reduces input value to [-PI, PI] range
double normalizeLon(double lon)
{
    while(lon >  M_PI) lon -= M_PI*2;
    while(lon < -M_PI) lon += M_PI*2;
    return lon;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    if(argc != 6 && argc != 7)
    {
        std::cerr << "Usage: " << argv[0] << " sphereRadiusInKM kmPerUnitValue inputFile outputFileName normalMapHeightInPixels [supersamplingLevel]\n";
        return 1;
    }

    const auto sphereRadius = std::stod(argv[1]);
    const auto kmPerUnit = std::stod(argv[2]);
    const char*const inFileName  = argv[3];
    const char*const outFileName = argv[4];

    const auto normalMapHeight = std::stoul(argv[5]);
    const auto normalMapWidth = 2*normalMapHeight;

    const unsigned numSamples = argc<7 ? 1 : std::max(1ul, std::stoul(argv[6]));

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

    const auto inData = reinterpret_cast<uint16_t*>(in.bits());
    const auto width  = in.width();
    const auto height = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof inData[0])
    {
        std::cerr << "Row stride of " << strideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto rowStride = strideInBytes / sizeof inData[0];

    std::vector<glm::dvec2> samplesCenter(numSamples);
    std::vector<glm::dvec2> samplesEast  (numSamples);
    std::vector<glm::dvec2> samplesNorth (numSamples);
    if(numSamples<=1)
        samplesCenter = samplesEast = samplesNorth = {glm::dvec2(0,0)};
    else
    {
        std::mt19937 mt(std::random_device{}());
        std::uniform_real_distribution<double> dist(-0.5,0.5);
        static const auto rng=[&](){return dist(mt);};

        for(auto& sample : samplesCenter)
            sample = glm::dvec2(rng(), rng());
        for(auto& sample : samplesEast)
            sample = glm::dvec2(rng(), rng());
        for(auto& sample : samplesNorth)
            sample = glm::dvec2(rng(), rng());
    }

    const size_t outBytesPerPixel = 3;
    using OutType = uchar;
    std::vector<OutType> outData(normalMapWidth*normalMapHeight*outBytesPerPixel);
    std::atomic<unsigned> rowsDone{0};
    const auto startTime = std::chrono::steady_clock::now();
    std::atomic_int numThreadsReportedFirstProgress{0};
    auto work = [inData, outData=outData.data(), numSamples, outBytesPerPixel,
                 normalMapWidth, normalMapHeight, width, height, rowStride,
                 kmPerUnit, sphereRadius, startTime,
                 &samplesCenter=std::as_const(samplesCenter),
                 &samplesEast=std::as_const(samplesEast),
                 &samplesNorth=std::as_const(samplesNorth),
                 &rowsDone, &numThreadsReportedFirstProgress]
                 (const size_t jMin, const size_t jMax)
    {
        auto time0 = std::chrono::steady_clock::now();
        size_t rowsDoneInThisThreadAfterLastUpdate = 0;
        for(size_t j = jMin; j < jMax; ++j)
        {
            const auto v = (0.5+(normalMapHeight - 1 - j)) / normalMapHeight;
            for(size_t i = 0; i < normalMapWidth; ++i)
            {
                const auto u = (0.5 + i) / normalMapWidth;
                const auto pixelPosInData = (i + j*normalMapWidth)*outBytesPerPixel;

                using namespace glm;
                using namespace std;

                double centerRadius = 0;
                double eastRadius = 0;
                double northRadius = 0;

                const double centerLon = (2*u-1)*M_PI;
                const double centerLat  = (2*v-1)*(M_PI/2);

                const double deltaLon = (2*M_PI/normalMapWidth) / cos(centerLat);
                const double eastLon = centerLon + deltaLon;
                const double eastLat = centerLat;

                const auto deltaLat = M_PI/normalMapHeight;
                const double northLat = centerLat + deltaLat;
                const double northLon = centerLon;

                for(unsigned sampleN = 0; sampleN < numSamples; ++sampleN)
                {
                    const double centerHeight = kmPerUnit*sample(inData, width, height, rowStride,
                                                                 normalizeLon(centerLon + samplesCenter[sampleN].x*deltaLon),
                                                                 centerLat + samplesCenter[sampleN].y*deltaLat);
                    centerRadius += sphereRadius + centerHeight;

                    const double eastHeight = kmPerUnit*sample(inData, width, height, rowStride,
                                                               normalizeLon(eastLon + samplesEast[sampleN].x*deltaLon),
                                                               eastLat + samplesEast[sampleN].y*deltaLat);
                    eastRadius += sphereRadius + eastHeight;

                    const double northHeight = kmPerUnit*sample(inData, width, height, rowStride,
                                                                normalizeLon(northLon + samplesNorth[sampleN].x*deltaLon),
                                                                northLat + samplesNorth[sampleN].y*deltaLat);
                    northRadius += sphereRadius + northHeight;
                }
                centerRadius /= numSamples;
                eastRadius   /= numSamples;
                northRadius  /= numSamples;

                const dvec3 centerDir = dvec3(cos(centerLon) * cos(centerLat),
                                              sin(centerLon) * cos(centerLat),
                                              sin(centerLat));
                const dvec3 centerPoint = centerRadius * centerDir;

                const dvec3 eastPoint = eastRadius * dvec3(cos(eastLon) * cos(eastLat),
                                                           sin(eastLon) * cos(eastLat),
                                                           sin(eastLat));
                const dvec3 northPoint = northRadius * dvec3(cos(northLon) * cos(northLat),
                                                             sin(northLon) * cos(northLat),
                                                             sin(northLat));
                const dvec3 deltaEast = eastPoint - centerPoint;
                const dvec3 deltaNorth = northPoint - centerPoint;

                const dvec3 normal = normalize(cross(deltaEast, deltaNorth));

                const dvec3 axis1 = normalize(cross(dvec3(0,0,1), centerDir));
                const dvec3 axis2 = normalize(cross(centerDir, axis1));
                const dvec3 axis3 = centerDir;

                const double normalA = dot(normal, axis1);
                const double normalB = dot(normal, axis2);
                const double normalC = dot(normal, axis3);

                constexpr auto max = std::numeric_limits<OutType>::max();
                outData[pixelPosInData + 0] = std::clamp(normalA+0.5,0.,1.) * max;
                outData[pixelPosInData + 1] = std::clamp(normalB+0.5,0.,1.) * max;
                outData[pixelPosInData + 2] = std::clamp(normalC    ,0.,1.) * max;
            }

            ++rowsDoneInThisThreadAfterLastUpdate;
            using namespace std::chrono;
            auto time1 = steady_clock::now();
            if(time1 - time0 > seconds(5))
            {
                // ETA is reported only by the first thread in the iteration.
                // This seems to provide the most accurate estimate. Once the
                // reporting thread is determined, it will report after each
                // iteration.
                thread_local bool isTimeReportingThread = false;
                thread_local bool reportedProgressFirstTime = false;
                if(!reportedProgressFirstTime)
                {
                    if(numThreadsReportedFirstProgress.fetch_add(1) == 0)
                        isTimeReportingThread = true;
                }

                rowsDone += rowsDoneInThisThreadAfterLastUpdate;
                rowsDoneInThisThreadAfterLastUpdate = 0;
                const auto progress = std::round(rowsDone * 10000. / normalMapHeight)/10000.;
                const auto usecElapsed = duration_cast<microseconds>(time1 - startTime).count();
                const long secToEnd = std::lround((1. - progress) / progress * usecElapsed * 1e-6);
                std::ostringstream ss;
                ss << progress*100 << "% done";
                // First ETA estimate is too far from reality, likely due
                // to overhead of thread startup, so skip first measurement.
                if(isTimeReportingThread && reportedProgressFirstTime)
                {
                    if(secToEnd > 60)
                    {
                        const auto min = secToEnd/60;
                        const auto sec = secToEnd - min*60;
                        ss << ", ETA: " << min << 'm' << sec << "s\n";
                    }
                    else
                    {
                        ss << ", ETA: " << secToEnd << "s\n";
                    }
                }
                else
                {
                    ss << '\n';
                }
                reportedProgressFirstTime = true;
                std::cerr << ss.str();
                time0 = time1;
            }
        }
    };

    auto time0 = std::chrono::steady_clock::now();
    const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for(size_t n = 0; n < numThreads; ++n)
    {
        const size_t jMin = normalMapHeight / numThreads * n;
        const size_t jMax = n+1 < numThreads ? normalMapHeight / numThreads * (n+1) : normalMapHeight;
        threads.emplace_back(work, jMin, jMax);
    }
    for(auto& thread : threads)
        thread.join();
    auto time1 = std::chrono::steady_clock::now();
    std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";

    std::cerr << "Saving image to " << outFileName << "... ";
    const QImage out(outData.data(), normalMapWidth, normalMapHeight,
                     normalMapWidth*outBytesPerPixel, QImage::Format_RGB888);
    if(!out.save(outFileName))
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
