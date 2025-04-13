#include <cmath>
#include <random>
#include <atomic>
#include <thread>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <QDir>
#include <QFile>
#include <QImage>
#include <QPainter>
#include <QRegularExpression>
#include "hips.hpp"
#include "timing.hpp"
#include "healpix.hpp"

namespace
{
constexpr double SPHERE_RADIUS = 1'737'400;

enum SectorOffset
{
    SO_NORTH,
    SO_SOUTH,
    SO_WEST,
    SO_LON0,
    SO_EAST,
    SO_LON180,
};

struct Tile
{
    QString sector;
    static constexpr int MAX_ALT_TILES_PER_SIDE = 2;
    double maxAltitude; // for the whole sector
    double maxAltitudes[MAX_ALT_TILES_PER_SIDE][MAX_ALT_TILES_PER_SIDE]; // per tile
    std::unique_ptr<QFile> file;
    const int16_t* data;
    ssize_t width, height;
    SectorOffset sectorOffset;
};

const char*const sectorNames[] = {
    "north",
    "south",
    "west",
    "lon0",
    "east",
    "lon180",
};

template<typename T> auto sqr(T const& x) { return x*x; }

//! Reduces input value to [-PI, PI] range
double normalizeLon(double lon)
{
    while(lon >  M_PI) lon -= M_PI*2;
    while(lon < -M_PI) lon += M_PI*2;
    return lon;
}

double normalizeLat(double lat)
{
    while(lat >  M_PI/2) lat -= M_PI*2;
    while(lat < -M_PI/2) lat += M_PI*2;
    return lat;
}

double fetch(Tile const& tile, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = std::clamp(requestedX, ssize_t(0), tile.width-1);
    const auto y = std::clamp(requestedY, ssize_t(0), tile.height-1);
    return tile.data[y*tile.width + x];
}

double sampleRect(Tile const& tile, const double u, const double v)
{
    const auto x = u * (tile.width - 1);
    const auto y = v * (tile.height - 1);

    const auto floorX = std::floor(x);
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetch(tile, floorX  , floorY);
    const auto pTopRight    = fetch(tile, floorX+1, floorY);
    const auto pBottomLeft  = fetch(tile, floorX  , floorY+1);
    const auto pBottomRight = fetch(tile, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;
    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

std::tuple<const Tile*,double,double> getTileAndUV(std::vector<Tile> const& tiles,
                                                   const glm::dvec3& dir)
{
    const auto x = dir[0];
    const auto y = dir[1];
    const auto z = dir[2];
    const auto ax = std::abs(x);
    const auto ay = std::abs(y);
    const auto az = std::abs(z);
    double u, v;

    const Tile* tile;
    if(ax >= ay && ax >= az)
    {
        // ±x is the dominant direction
        if(x > 0)
        {
            // {1, u*2-1, v*2-1} => lon = 0
            tile = &tiles[SO_LON0];
            u = y / x * 0.5 + 0.5;
            v = z / x * 0.5 + 0.5;
        }
        else
        {
            // {-1, -(u*2-1), v*2-1} => lon = 180°
            tile = &tiles[SO_LON180];
            u = -y / ax * 0.5 + 0.5;
            v =  z / ax * 0.5 + 0.5;
        }
    }
    else if(ay >= ax && ay >= az)
    {
        // ±y is the dominant direction
        if(y > 0)
        {
            // {-(u*2-1), 1, v*2-1} => east
            tile = &tiles[SO_EAST];
            u = -x / y * 0.5 + 0.5;
            v =  z / y * 0.5 + 0.5;
        }
        else
        {
            // {u*2-1, -1, v*2-1} => west
            tile = &tiles[SO_WEST];
            u = x / ay * 0.5 + 0.5;
            v = z / ay * 0.5 + 0.5;
        }
    }
    else
    {
        // ±z is the dominant direction
        if(z > 0)
        {
            // {u*2-1, v*2-1, 1} => north
            tile = &tiles[SO_NORTH];
            u = x / z * 0.5 + 0.5;
            v = y / z * 0.5 + 0.5;
        }
        else
        {
            // {u*2-1, -(v*2-1), -1} => south
            tile = &tiles[SO_SOUTH];
            u =  x / az * 0.5 + 0.5;
            v = -y / az * 0.5 + 0.5;
        }
    }
    return {tile, u, v};
}

double sample(std::vector<Tile> const& tiles, const glm::dvec3& dir)
{
    const auto [tile, u, v] = getTileAndUV(tiles, dir);
    return sampleRect(*tile, u, v);
}

struct SubsectorId
{
    unsigned sectorOffset;
    unsigned subsectorI, subsectorJ;
    double u, v;
    bool operator==(SubsectorId const& r) const
    {
        return sectorOffset == r.sectorOffset &&
               subsectorI == r.subsectorI &&
               subsectorJ == r.subsectorJ;
    }
};

SubsectorId dirToSubsector(const std::vector<Tile>& tiles, glm::dvec3 const& dir)
{
    const auto [tile, u, v] = getTileAndUV(tiles, dir);
    const unsigned x = std::lround(u * (tile->width - 1));
    const unsigned y = std::lround(v * (tile->height - 1));

    const unsigned subsecI = x * Tile::MAX_ALT_TILES_PER_SIDE / tile->width;
    const unsigned subsecJ = y * Tile::MAX_ALT_TILES_PER_SIDE / tile->height;
    return {tile->sectorOffset, subsecI, subsecJ, u, v};
}

double subsectorsMaxRadiusSquared(const std::vector<Tile>& heightMapTiles,
                                  const SubsectorId srcSector,
                                  const SubsectorId finalSector,
                                  const double globalMax)
{
    if(srcSector.sectorOffset == finalSector.sectorOffset &&
       srcSector.sectorOffset != SO_NORTH && srcSector.sectorOffset != SO_SOUTH)
    {
        // In the simple (non-equiangular) cubemap projection all great circles of the concentric
        // sphere project to straight lines on each face. So we can assume that vertical and
        // horizontal line segments that start in different tiles of the same equatorial face only
        // cover these two tiles (the additional assumption being that these segments are shorter
        // than the tile widths). Our north-south and east-west lines are respectively vertical and
        // horizontal.
        const auto sec = srcSector.sectorOffset;
        const auto i1 = srcSector.subsectorI, j1 = srcSector.subsectorJ;
        const auto i2 = finalSector.subsectorI, j2 = finalSector.subsectorJ;
        const auto maxAlt = std::max(heightMapTiles[sec].maxAltitudes[i1][j1],
                                     heightMapTiles[sec].maxAltitudes[i2][j2]);
        return sqr(SPHERE_RADIUS + maxAlt);
    }

    // Even if we haven't determined the exact set of tiles used, we
    // can optimize for the case of the same face being sampled.
    if(srcSector.sectorOffset == finalSector.sectorOffset)
    {
        const auto maxAlt = heightMapTiles[srcSector.sectorOffset].maxAltitude;
        return sqr(SPHERE_RADIUS + maxAlt);
    }

    return globalMax;
}

Tile readTile(QString const& inDir, QString const& sector)
{
    Tile tile;
    tile.sector = sector;
    const auto filename = inDir + "/cube-face-" + sector + ".dat";
    tile.file.reset(new QFile(filename));
    if(!tile.file->open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << filename.toStdString()
                  << " for reading: "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    uint16_t formatVersion;
    if(tile.file->read(reinterpret_cast<char*>(&formatVersion), sizeof formatVersion) != sizeof formatVersion)
    {
        std::cerr << "Failed to read cube map format version "
                  << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    constexpr uint16_t SUPPORTED_VERSION = 2;
    if(formatVersion != SUPPORTED_VERSION)
    {
        std::cerr << filename.toStdString() << ": cube map format version " << formatVersion
                  << " differs from the expected version " << SUPPORTED_VERSION << ".\n";
        return {};
    }
    ssize_t cubeMapSide = -1;
    if(tile.file->read(reinterpret_cast<char*>(&cubeMapSide), sizeof cubeMapSide) != sizeof cubeMapSide)
    {
        std::cerr << "Failed to read cube map side length of "
                  << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    int16_t maxAltitudes[Tile::MAX_ALT_TILES_PER_SIDE][Tile::MAX_ALT_TILES_PER_SIDE];
    if(tile.file->read(reinterpret_cast<char*>(&maxAltitudes), sizeof maxAltitudes) != sizeof maxAltitudes)
    {
        std::cerr << "Failed to read cube map maximum value "
                  << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    std::cerr << "Maximum altitudes in the sector:";
    tile.maxAltitude = -INFINITY;
    for(int j = 0; j < Tile::MAX_ALT_TILES_PER_SIDE; ++j)
    {
        for(int i = 0; i < Tile::MAX_ALT_TILES_PER_SIDE; ++i)
        {
            tile.maxAltitudes[i][j] = maxAltitudes[i][j];
            if(tile.maxAltitudes[i][j] > tile.maxAltitude)
                tile.maxAltitude = tile.maxAltitudes[i][j];
            std::cerr << (i==0 && j==0 ? " " : ", ") << tile.maxAltitudes[i][j];
        }
    }
    std::cerr << "\n";
    std::cerr << "Maximum altitude in the whole sector: " << tile.maxAltitude << " m\n";

    tile.width = cubeMapSide;
    tile.height= cubeMapSide;
    if(tile.width==0 || tile.height==0)
    {
        std::cerr << "Failed to find image dimensions\n";
        return {};
    }
    std::cerr << "Image dimensions for sector " << sector.toStdString()
              << ": " << tile.width << u8"×" << tile.height << "\n";

    const qint64 sizeToRead = tile.width * tile.height * sizeof tile.data[0];
    const qint64 headerSize = tile.file->pos();
    tile.data = reinterpret_cast<const int16_t*>(tile.file->map(headerSize, sizeToRead));
    if(!tile.data)
    {
        std::cerr << "Failed to read " << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    return tile;
}

glm::dvec4 computeHorizons(const double raySourceLon, const double raySourceLat,
                           const std::vector<Tile>& heightMapTiles,
                           const double resolutionAtEquator,
                           const double maxRadiusSquaredOrig)
{
    using namespace glm;

    const double deltaLatLon = 2*M_PI/resolutionAtEquator;
    const double eastLon = raySourceLon + deltaLatLon;
    const double eastLat = raySourceLat;

    const double northLat = raySourceLat + deltaLatLon;
    const double northLon = raySourceLon;

    const dvec3 zenithAtRaySource = dvec3(cos(raySourceLon) * cos(raySourceLat),
                                          sin(raySourceLon) * cos(raySourceLat),
                                          sin(raySourceLat));
    const dvec3 zenithAtEast = dvec3(cos(eastLon) * cos(eastLat),
                                     sin(eastLon) * cos(eastLat),
                                     sin(eastLat));
    const dvec3 zenithAtNorth = dvec3(cos(northLon) * cos(northLat),
                                      sin(northLon) * cos(northLat),
                                      sin(northLat));
    const double raySourceHeight = sample(heightMapTiles, zenithAtRaySource);
    const double raySourceRadius = SPHERE_RADIUS + raySourceHeight;
    const dvec3 raySourcePoint = raySourceRadius * zenithAtRaySource;

    const double eastHeight = sample(heightMapTiles, zenithAtEast);
    const double eastRadius = SPHERE_RADIUS + eastHeight;
    const dvec3 eastPoint = eastRadius * zenithAtEast;

    const double northHeight = sample(heightMapTiles, zenithAtNorth);
    const double northRadius = SPHERE_RADIUS + northHeight;
    const dvec3 northPoint = northRadius * zenithAtNorth;

    const dvec3 deltaEast = eastPoint - raySourcePoint;
    const dvec3 deltaNorth = northPoint - raySourcePoint;

    glm::dvec4 horiz;
    unsigned rayIndex = 0;
    for(dvec3 rayDir : {deltaNorth, deltaEast, -deltaNorth, -deltaEast})
    {
        const auto srcSector = dirToSubsector(heightMapTiles, raySourcePoint);
        const auto maxRayLength = std::sqrt(maxRadiusSquaredOrig - sqr(raySourceRadius));
        constexpr double safetyMargin = 1.01; // to avoid missing a border nearby
        const auto finalRayPoint = raySourcePoint + normalize(rayDir) * (maxRayLength * safetyMargin);
        const auto finalSector = dirToSubsector(heightMapTiles, finalRayPoint);
        const auto maxRadiusSquared = subsectorsMaxRadiusSquared(heightMapTiles, srcSector, finalSector, maxRadiusSquaredOrig);

        double sinHorizonElevation = -M_PI/2;
        for(dvec3 rayPoint = raySourcePoint + rayDir; dot(rayPoint, rayPoint) <= maxRadiusSquared; rayPoint += rayDir)
        {
            const double altitude = sample(heightMapTiles, rayPoint);
            const dvec3 zenithAtRayPoint = normalize(rayPoint);
            const dvec3 pointAtAlt = zenithAtRayPoint*(SPHERE_RADIUS+altitude);
            const double sinElev = dot(normalize(pointAtAlt-raySourcePoint), zenithAtRaySource);
            if(sinElev > sinHorizonElevation)
            {
                sinHorizonElevation = sinElev;
                // Tilt the ray up to point above the current horizon estimate
                // rayDir = length(rayDir) * normalize(pointAtAlt - raySourcePoint);
                const dvec3 newRayUnnorm = pointAtAlt - raySourcePoint;
                rayDir = newRayUnnorm * std::sqrt(dot(rayDir, rayDir) / dot(newRayUnnorm, newRayUnnorm));
                rayPoint = pointAtAlt;
            }
        }
        horiz[rayIndex] = sinHorizonElevation;
        ++rayIndex;
    }
    return horiz;
}

void fillHorizonsFace(const int order, const int pix, const std::vector<Tile>& heightMapTiles,
                      const int channelsPerPixel, const double resolutionAtEquator,
                      const double maxRadiusSquared, uint8_t* outData)
{
    const unsigned nside = 1u << order;
    int ix, iy, face;
    healpix_nest2xyf(nside, pix, &ix, &iy, &face);

    for(int y = 0; y < HIPS_TILE_SIZE; ++y)
    {
        for(int x = 0; x < HIPS_TILE_SIZE; ++x)
        {
            double theta, phi;
            healpix_xyf2ang(nside * HIPS_TILE_SIZE,
                            ix * HIPS_TILE_SIZE + x, iy * HIPS_TILE_SIZE + y,
                            face, &theta, &phi);
            const double latitude = normalizeLat(M_PI/2 - theta);
            const double longitude = normalizeLon(M_PI+phi);
            assert(-M_PI <= longitude && longitude <= M_PI);
            assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

            // HiPS coordinates are swapped, because they were designed from the "look from
            // sphere center" perspective, while we are making a look from outside a planet.
            const int i = y, j = x;

            const int pixelPosInOutData = (i + j*HIPS_TILE_SIZE)*channelsPerPixel;
            const auto sinHorizElevs = computeHorizons(longitude, latitude, heightMapTiles,
                                                       resolutionAtEquator,
                                                       maxRadiusSquared);
            const auto v = sign(sinHorizElevs)*sqrt(abs(sinHorizElevs));
            constexpr auto max = std::numeric_limits<uint8_t>::max();
            for(int i = 0; i < channelsPerPixel; ++i)
                outData[pixelPosInOutData + i] = std::clamp(v[i]/2+0.5, 0.,1.) * max;
        }
    }
}

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} inDir outDir";
    s << R"(
Options:
 -h, --help                 This help message
 --format FORMAT            Set hips_tile_format to FORMAT (only one value is supported)
 --title TITLE              Set obs_title to TITLE
 -d, --desc DESCRIPTION     Set obs_description to DESCRIPTION
 --frame FRAME              Set hips_frame to FRAME
 --creator CREATOR          Set hips_creator to CREATOR
 --my-copyright COPYRIGHT   Set hips_copyright to COPYRIGHT
 --orig-copyright COPYRIGHT Set obs_copyright to COPYRIGHT
 -s, --status STATUS        Set hips_status to STATUS
)";
    return ret;
}

}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    QString inDir;
    QString outDir;
    QString imgFormat = "png";
    QString surveyTitle;
    QString description;
    QString frame;
    QString creator;
    QString hips_copyright;
    QString obs_copyright;
    QString hipsStatus = "public mirror clonable";

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

        if(arg == "--format")
        {
            GO_TO_PARAM();
            imgFormat = argv[n];
        }
        else if(arg == "--title")
        {
            GO_TO_PARAM();
            surveyTitle = argv[n];
        }
        else if(arg == "-d" || arg == "--desc" || arg == "--description")
        {
            GO_TO_PARAM();
            description = argv[n];
        }
        else if(arg == "--frame")
        {
            GO_TO_PARAM();
            frame = argv[n];
        }
        else if(arg == "--creator")
        {
            GO_TO_PARAM();
            creator = argv[n];
        }
        else if(arg == "--my-copyright")
        {
            GO_TO_PARAM();
            hips_copyright = argv[n];
        }
        else if(arg == "--orig-copyright")
        {
            GO_TO_PARAM();
            obs_copyright = argv[n];
        }
        else if(arg == "-s" || arg == "--status")
        {
            GO_TO_PARAM();
            hipsStatus = argv[n];
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(inDir.isEmpty())
    {
        std::cerr << "Input file not specified\n";
        return usage(argv[0], 1);
    }
    if(outDir.isEmpty())
    {
        std::cerr << "Output directory not specified\n";
        return usage(argv[0], 1);
    }

    QString finalExt;
    if(imgFormat == "jpeg")
        finalExt = "jpg";
    else if(imgFormat == "webp" || imgFormat == "bmp" || imgFormat == "png" || imgFormat == "tiff")
        finalExt = imgFormat;
    else
        throw std::runtime_error("Unexpected output image format: " + imgFormat.toStdString());

    std::vector<Tile> heightMapTiles;
    for(const QString& sector : sectorNames)
    {
        const auto sectorOffset = heightMapTiles.size();
        heightMapTiles.push_back(readTile(inDir, sector));
        heightMapTiles.back().sectorOffset = static_cast<SectorOffset>(sectorOffset);
    }
    for(unsigned i = 1; i < heightMapTiles.size(); ++i)
    {
        if(heightMapTiles[i].width != heightMapTiles[i-1].width || heightMapTiles[i].height != heightMapTiles[i-1].height)
        {
            throw std::runtime_error(QString(u8"With and height of sector %1 don't match that of %2: %3×%4 vs %5×%6")
                                        .arg(heightMapTiles[i].sector).arg(heightMapTiles[i-1].sector)
                                        .arg(heightMapTiles[i].width).arg(heightMapTiles[i].height)
                                        .arg(heightMapTiles[i-1].width).arg(heightMapTiles[i-1].height)
                                        .toStdString());
        }
    }

    double maxAltitude = -INFINITY;
    for(const auto& tile : heightMapTiles)
    {
        for(int j = 0; j < Tile::MAX_ALT_TILES_PER_SIDE; ++j)
        {
            for(int i = 0; i < Tile::MAX_ALT_TILES_PER_SIDE; ++i)
            {
                if(tile.maxAltitudes[i][j] > maxAltitude)
                    maxAltitude = tile.maxAltitudes[i][j];
            }
        }
    }
    const double maxRadiusSquared = sqr(SPHERE_RADIUS+maxAltitude);

    constexpr double resolutionAtEquator = 27291*4; // original data pixels per 360°

    if(!QDir().mkpath(outDir))
        throw std::runtime_error("Failed to create directory \""+outDir.toStdString()+'"');

    const int orderMax = std::ceil(std::log2(resolutionAtEquator / (4. * HIPS_TILE_SIZE * M_SQRT2)));

    hipsSaveProperties(outDir, orderMax, imgFormat, surveyTitle, "planet-horizon",
                       description, frame, obs_copyright, hips_copyright, creator, hipsStatus);

    // First create the tiles of the deepest level
    std::cerr << "Creating tiles of order " << orderMax << "...\n";
    {
        const int absolutePixMax = 12 * (1 << (2 * orderMax));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0},  itemsSkipped{0};
        const auto startTime = std::chrono::steady_clock::now();
        std::atomic_int currPix{0};
        auto work = [absolutePixMax,orderMax,outDir,startTime,&heightMapTiles = std::as_const(heightMapTiles),
                     &currPix,resolutionAtEquator,maxRadiusSquared,
                     &numThreadsReportedFirstProgress,&itemsDone,&itemsSkipped]()
        {
            constexpr int channelsPerPixel = 4;
            std::vector<uint8_t> data(HIPS_TILE_SIZE * HIPS_TILE_SIZE * channelsPerPixel);

            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = currPix.fetch_add(1, std::memory_order_relaxed);
                pix < absolutePixMax;
                pix = currPix.fetch_add(1, std::memory_order_relaxed))
            {
                const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(orderMax).arg((pix / 10000) * 10000);
                if(!QDir().mkpath(outPath))
                    throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
                const auto fileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(hipsInitialExt);
                if(QFileInfo(fileName).exists())
                {
                    ++itemsSkipped;
                    std::cerr << "Skipping existing "+fileName.mid(outDir.size()+1).toStdString()+"\n";
                    continue;
                }
                    fillHorizonsFace(orderMax, pix, heightMapTiles, channelsPerPixel, resolutionAtEquator,
                                     maxRadiusSquared, data.data());
                const QImage out(data.data(), HIPS_TILE_SIZE, HIPS_TILE_SIZE, channelsPerPixel*HIPS_TILE_SIZE,
                                 QImage::Format_RGBA8888);
                if(!out.save(fileName, nullptr, 100))
                    throw std::runtime_error("Failed to save output file " + fileName.toStdString());

                handleProgressReporting(absolutePixMax - itemsSkipped, startTime, time0, numThreadsReportedFirstProgress,
                                        itemsDoneInThisThreadAfterLastUpdate, itemsDone);
            }
        };
        const auto time0 = std::chrono::steady_clock::now();
        const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        for(size_t n = 0; n < numThreads; ++n)
            threads.emplace_back(work);
        for(auto& thread : threads)
            thread.join();
        auto time1 = std::chrono::steady_clock::now();
        std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";
    }

    generateLowerOrderTiles(orderMax, outDir);
    convertTiles(finalExt, imgFormat, orderMax, outDir);
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
