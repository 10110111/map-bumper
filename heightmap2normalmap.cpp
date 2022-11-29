#include <cmath>
#include <chrono>
#include <random>
#include <iostream>
#include <algorithm>

#include <glm/glm.hpp>

#include <QImage>

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

    const auto data = reinterpret_cast<uint16_t*>(in.bits());
    const auto width  = in.width();
    const auto height = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof data[0])
    {
        std::cerr << "Row stride of " << strideInBytes << " bytes is not a multiple of a pixel, this is not supported\n";
        return 1;
    }
    const auto rowStride = strideInBytes / sizeof data[0];

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
    std::vector<double> outData(normalMapWidth*normalMapHeight*outBytesPerPixel);
    auto time0 = std::chrono::steady_clock::now();
    for(size_t i = 0; i < normalMapWidth; ++i)
    {
        const auto u = (0.5 + i) / normalMapWidth;
        for(size_t j = 0; j < normalMapHeight; ++j)
        {
            const auto v = (0.5+(normalMapHeight - 1 - j)) / normalMapHeight;
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
                const double centerHeight = kmPerUnit*sample(data, width, height, rowStride,
                                                             normalizeLon(centerLon + samplesCenter[sampleN].x*deltaLon),
                                                                          centerLat + samplesCenter[sampleN].y*deltaLat);
                centerRadius += sphereRadius + centerHeight;

                const double eastHeight = kmPerUnit*sample(data, width, height, rowStride,
                                                           normalizeLon(eastLon + samplesEast[sampleN].x*deltaLon),
                                                                        eastLat + samplesEast[sampleN].y*deltaLat);
                eastRadius += sphereRadius + eastHeight;

                const double northHeight = kmPerUnit*sample(data, width, height, rowStride,
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

            outData[pixelPosInData + 0] = normalA+0.5;
            outData[pixelPosInData + 1] = normalB+0.5;
            outData[pixelPosInData + 2] = normalC;
        }
        auto time1 = std::chrono::steady_clock::now();
        if(time1 - time0 > std::chrono::seconds(5))
        {
            std::cerr << std::round((i+1) * 10000. / normalMapWidth)/100. << "% done\n";
            time0 = time1;
        }
    }

    using OutType = uchar;
    std::vector<OutType> outBits;
    outBits.reserve(outData.size());
    for(auto v : outData)
        outBits.push_back(std::clamp(v,0.,1.)*std::numeric_limits<OutType>::max());

    const QImage out(outBits.data(), normalMapWidth, normalMapHeight, normalMapWidth*outBytesPerPixel, QImage::Format_RGB888);
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
