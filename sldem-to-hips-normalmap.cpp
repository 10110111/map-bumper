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
#include <proj.h>
#include <geotiff/xtiffio.h>
#include <geotiff/geotiff.h>
#include <geotiff/geovalues.h>
#include <geotiff/geo_normalize.h>
#include "hips.hpp"
#include "timing.hpp"
#include "healpix.hpp"
#include "cubemap.hpp"

namespace
{

constexpr double southLDEMmaxLat = -60 * (M_PI/180);

struct TIFFTile
{
    std::unique_ptr<float[]> data;
    ssize_t width, height;
    GTIF* gtiff;
    GTIFDefn defn;
    PJ* pj = nullptr;
    PJ_CONTEXT* ctx = nullptr;
    ~TIFFTile()
    {
        if(pj) proj_destroy(pj);
        if(ctx) proj_context_destroy(ctx);
    }
};

struct Tile
{
    QString sector;
    double sphereRadius;
    double metersPerUnit;
    std::unique_ptr<QFile> file;
    const float* data;
    ssize_t width, height;
};

enum SectorOffset
{
    // LDEM_512
    L_N45_180, L_N45_270, L_N45_000, L_N45_090,
    L_N00_180, L_N00_270, L_N00_000, L_N00_090,
    L_S45_180, L_S45_270, L_S45_000, L_S45_090,
    L_S90_180, L_S90_270, L_S90_000, L_S90_090,
    // SLDEM2015 512ppd
    SL_N30_180, SL_N30_225, SL_N30_270, SL_N30_315, SL_N30_000, SL_N30_045, SL_N30_090, SL_N30_135,
    SL_N00_180, SL_N00_225, SL_N00_270, SL_N00_315, SL_N00_000, SL_N00_045, SL_N00_090, SL_N00_135,
    SL_S30_180, SL_S30_225, SL_S30_270, SL_S30_315, SL_S30_000, SL_S30_045, SL_S30_090, SL_S30_135,
    SL_S60_180, SL_S60_225, SL_S60_270, SL_S60_315, SL_S60_000, SL_S60_045, SL_S60_090, SL_S60_135,

    LDEM_HORIZ_TILES_COUNT = L_N00_180 - L_N45_180,
    SLDEM_HORIZ_TILES_COUNT = SL_N00_180 - SL_N30_180,
    VERT_TILES_COUNT = 4,
    LDEM_FIRST = L_N45_180,
    LDEM_LAST = L_S90_090,
    SLDEM_FIRST = SL_N30_180,
    SLDEM_LAST = SL_S60_135,
    FIRST_SECTOR = LDEM_FIRST,
    LAST_SECTOR = SL_S60_135,
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

double fetch(std::vector<Tile> const& data, const ssize_t requestedX, const ssize_t requestedY, const ssize_t resolutionAtEquator, const bool useSLDEM)
{
    const auto totalWidth = resolutionAtEquator;
    const auto totalHeight = totalWidth / 2;
    const auto x = (requestedX+totalWidth) % totalWidth;
    auto y = std::clamp(requestedY, ssize_t(0), totalHeight-1);

    if(useSLDEM)
    {
        const auto upperBorderPos = (data[L_N45_180].height - data[SL_N30_180].height) * (VERT_TILES_COUNT / 2);
        const auto lowerBorderPos = data[SLDEM_FIRST].height * VERT_TILES_COUNT + upperBorderPos;
        // If the data point is not represented in SLDEM2015, use LDEM
        if(y < upperBorderPos)
        {
            const auto tileWidth = data[LDEM_FIRST].width;
            const auto xOffsetIndex = x / tileWidth;
            return data[LDEM_FIRST + xOffsetIndex].data[y * tileWidth + x % tileWidth];
        }
        if(y >= lowerBorderPos)
        {
            const auto tileWidth = data[L_S90_180].width;
            const auto xOffsetIndex = x / tileWidth;
            const auto y0 = data[L_S90_180].height * 3;
            return data[L_S90_180 + xOffsetIndex].data[(y - y0) * tileWidth + x % tileWidth];
        }
        // So the point is represented in SLDEM2015, fetch it
        const auto tileWidth = data[SLDEM_FIRST].width;
        const auto tileHeight= data[SLDEM_FIRST].height;
        const auto xOffsetIndex = x / tileWidth;
        const auto yInSLDEM = y - upperBorderPos;
        const auto yOffsetIndex = yInSLDEM / tileHeight;
        const auto yInTile = yInSLDEM % tileHeight;
        const auto& tile = data[SLDEM_FIRST + yOffsetIndex * SLDEM_HORIZ_TILES_COUNT + xOffsetIndex];
        return tile.data[yInTile * tileWidth + x % tileWidth];
    }
    else
    {
        const auto tileWidth = data[LDEM_FIRST].width;
        const auto tileHeight= data[LDEM_FIRST].height;
        const auto xOffsetIndex = x / tileWidth;
        const auto yOffsetIndex = y / tileHeight;
        const auto yInTile = y % tileHeight;
        const auto& tile = data[LDEM_FIRST + yOffsetIndex * LDEM_HORIZ_TILES_COUNT + xOffsetIndex];
        return tile.data[yInTile * tileWidth + x % tileWidth];
    }
}

double fetchFromRect(float const*const data, const ssize_t width, const ssize_t height,
                     const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = std::clamp(requestedX, ssize_t(0), width-1);
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);
    return data[y*width + x];
}

double sampleSouthLDEM(TIFFTile const& southLDEM, const double longitude, const double latitude)
{
    PJ_COORD coord;
    coord.xyzt.x = longitude * (180/M_PI);
    coord.xyzt.y = latitude * (180/M_PI);
    coord.xyzt.z = 0;
    coord.xyzt.t = 0;
    coord = proj_trans(southLDEM.pj, PJ_FWD, coord);
    double x = coord.xyzt.x;
    double y = coord.xyzt.y;

    GTIFPCSToImage(southLDEM.gtiff, &x, &y);
    if(x < 0 || y < 0 || x >= southLDEM.width || y >= southLDEM.height)
    {
        // Out of image, must use other data for this location
        return NAN;
    }

    const auto floorX = std::floor(x);
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetchFromRect(southLDEM.data.get(), southLDEM.width, southLDEM.height, floorX  , floorY);
    const auto pTopRight    = fetchFromRect(southLDEM.data.get(), southLDEM.width, southLDEM.height, floorX+1, floorY);
    const auto pBottomLeft  = fetchFromRect(southLDEM.data.get(), southLDEM.width, southLDEM.height, floorX  , floorY+1);
    const auto pBottomRight = fetchFromRect(southLDEM.data.get(), southLDEM.width, southLDEM.height, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;
    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

double sample(std::vector<Tile> const& data, TIFFTile const& southLDEM, const double resolutionAtEquator, double longitude, double latitude, const bool useSLDEM)
{
    if(latitude > M_PI/2 || latitude < -M_PI/2)
    {
        latitude = -latitude;
        longitude += M_PI;
    }

    if(latitude < southLDEMmaxLat)
    {
        const auto value = sampleSouthLDEM(southLDEM, longitude, latitude);
        const auto unitsPerKM = 1e-3;
        if(!std::isnan(value)) return value * unitsPerKM;
    }

    const auto deltaLon = 2.*M_PI/resolutionAtEquator;
    const auto firstLon = (1.-resolutionAtEquator)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -M_PI/(resolutionAtEquator/2-1);
    const auto firstLat = (1.-resolutionAtEquator/2)/2. * deltaLat;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetch(data, floorX  , floorY  , resolutionAtEquator, useSLDEM);
    const auto pTopRight    = fetch(data, floorX+1, floorY  , resolutionAtEquator, useSLDEM);
    const auto pBottomLeft  = fetch(data, floorX  , floorY+1, resolutionAtEquator, useSLDEM);
    const auto pBottomRight = fetch(data, floorX+1, floorY+1, resolutionAtEquator, useSLDEM);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

Tile readTile(QString const& inDir, QString const& sector)
{
    Tile tile;
    tile.sector = sector;
    const auto filenameBase = inDir + "/" + sector + "_FLOAT";

    const auto labelFileName = filenameBase + ".LBL";
    QFile file(labelFileName);
    if(!file.open(QFile::ReadOnly))
    {
        throw std::runtime_error(QString("Failed to open %1 for reading: %2")
                                    .arg(labelFileName).arg(file.errorString()).toStdString());
    }
    const auto label = file.readAll();
    if(label.isEmpty())
    {
        throw std::runtime_error(QString("Failed to read label file %1: %2")
                                    .arg(labelFileName).arg(file.errorString()).toStdString());
    }

    const auto filename = filenameBase + ".IMG";
    tile.file.reset(new QFile(filename));
    if(!tile.file->open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << filename.toStdString()
                  << " for reading: "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }

    tile.width = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();
    tile.height= QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();
    if(tile.width==0 || tile.height==0)
    {
        std::cerr << "Failed to find image dimensions\n";
        return {};
    }
    std::cerr << "Image dimensions for sector " << sector.toStdString()
              << ": " << tile.width << u8"×" << tile.height << "\n";
    const auto radiusA = QRegularExpression(R"(\n\s*A_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();
    const auto radiusB = QRegularExpression(R"(\n\s*B_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();
    const auto radiusC = QRegularExpression(R"(\n\s*C_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();
    if(radiusA <= 0 || radiusA != radiusB || radiusB != radiusC)
        throw std::runtime_error(QString("Bad reference sphere radii found. Expected three equal positive radii, but found %1, %2, %3")
                                    .arg(radiusA).arg(radiusB).arg(radiusC).toStdString());
    tile.sphereRadius = 1000 * radiusA;

    const auto sampleOffset = QRegularExpression(R"(\n\s*\bOFFSET\s*=\s*([0-9.]+)\s*\n)").match(label).captured(1).toDouble();
    if(std::abs(sampleOffset - radiusA) > 1e-6)
        throw std::runtime_error(QString("Unexpected sample offset %1 (expected %2)").arg(sampleOffset).arg(radiusA).toStdString());

    const auto sampleType = QRegularExpression(R"(\n\s*\bSAMPLE_TYPE\s*=\s*([^\s]+)\s*\n)").match(label).captured(1);
    const auto sampleBits = QRegularExpression(R"(\n\s*\bSAMPLE_BITS\s*=\s*([^\s]+)\s*\n)").match(label).captured(1);
    if(sampleType != "PC_REAL" || sampleBits != "32")
    {
        throw std::runtime_error(QString("Sample type %1 with %2 bits per sample isn't supported. Only (little-endian) float32 is.")
                                    .arg(sampleType).arg(sampleBits).toStdString());
    }

    const auto sampleUnit = QRegularExpression(R"(\n\s*\bUNIT\s*=\s*([^\s]+)\b\s*\n)").match(label).captured(1);
    if(sampleUnit != "KILOMETER")
        throw std::runtime_error(("Unexpected sample unit: "+sampleUnit).toStdString());
    const auto sampleScalingFactor = QRegularExpression(R"(\bSCALING_FACTOR\s*=\s*([0-9.]+)\s*\n)").match(label).captured(1).toDouble();
    tile.metersPerUnit = sampleScalingFactor * 1000;

    const auto projType = QRegularExpression(R"re(\n\s*MAP_PROJECTION_TYPE\s*=\s*"([^"]+)"\s*\n)re").match(label).captured(1);
    if(projType != "SIMPLE CYLINDRICAL")
        throw std::runtime_error(QString("Unexpected projection type: %s, while only simple cylindrical is supported.").arg(projType).toStdString());

    const qint64 sizeToRead = tile.width * tile.height * sizeof tile.data[0];
    tile.data = reinterpret_cast<const float*>(tile.file->map(0, sizeToRead));
    if(!tile.data)
    {
        std::cerr << "Failed to read " << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    return tile;
}

glm::dvec3 computeNormal(const double centerLon, const double centerLat,
                         const std::vector<Tile>& heightMapTiles, const TIFFTile& southLDEM,
                         const double resolutionAtEquator,
                         const double metersPerUnit, const double sphereRadius, const bool useSLDEM)
{
    using namespace glm;

    const double deltaLon = (2*M_PI/resolutionAtEquator) / cos(centerLat);
    const double eastLon = centerLon + deltaLon;
    const double eastLat = centerLat;
    const double westLon = centerLon - deltaLon;
    const double westLat = centerLat;

    const double deltaLat = M_PI/resolutionAtEquator;
    const double northLat = centerLat + deltaLat;
    const double northLon = centerLon;
    const double southLat = centerLat - deltaLat;
    const double southLon = centerLon;

    const double eastHeight = metersPerUnit*sample(heightMapTiles, southLDEM, resolutionAtEquator, eastLon, eastLat, useSLDEM);
    const double eastRadius = sphereRadius + eastHeight;

    const double westHeight = metersPerUnit*sample(heightMapTiles, southLDEM, resolutionAtEquator, westLon, westLat, useSLDEM);
    const double westRadius = sphereRadius + westHeight;

    const double northHeight = metersPerUnit*sample(heightMapTiles, southLDEM, resolutionAtEquator, northLon, northLat, useSLDEM);
    const double northRadius = sphereRadius + northHeight;

    const double southHeight = metersPerUnit*sample(heightMapTiles, southLDEM, resolutionAtEquator, southLon, southLat, useSLDEM);
    const double southRadius = sphereRadius + southHeight;

    const dvec3 centerDir = dvec3(cos(centerLon) * cos(centerLat),
                                  sin(centerLon) * cos(centerLat),
                                  sin(centerLat));

    const dvec3 eastPoint = eastRadius * dvec3(cos(eastLon) * cos(eastLat),
                                               sin(eastLon) * cos(eastLat),
                                               sin(eastLat));
    const dvec3 westPoint = westRadius * dvec3(cos(westLon) * cos(westLat),
                                               sin(westLon) * cos(westLat),
                                               sin(westLat));
    const dvec3 northPoint = northRadius * dvec3(cos(northLon) * cos(northLat),
                                                 sin(northLon) * cos(northLat),
                                                 sin(northLat));
    const dvec3 southPoint = southRadius * dvec3(cos(southLon) * cos(southLat),
                                                 sin(southLon) * cos(southLat),
                                                 sin(southLat));
    const dvec3 deltaEast = eastPoint - westPoint;
    const dvec3 deltaNorth = northPoint - southPoint;

    const dvec3 normal = normalize(cross(deltaEast, deltaNorth));

    const dvec3 axis1 = normalize(cross(dvec3(0,0,1), centerDir));
    const dvec3 axis2 = normalize(cross(centerDir, axis1));
    const dvec3 axis3 = centerDir;

    const double normalA = dot(normal, axis1);
    const double normalB = dot(normal, axis2);
    const double normalC = dot(normal, axis3);
    return {normalA, normalB, normalC};
}

void fillFace(const int order, const int pix, const std::vector<Tile>& heightMapTiles,
              const TIFFTile& southLDEM,
              const int channelsPerPixel, const double resolutionAtEquator,
              const double metersPerUnit, const double sphereRadius, const bool useSLDEM,
              uint8_t* outData)
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
            const auto normal = computeNormal(longitude, latitude, heightMapTiles, southLDEM,
                                              resolutionAtEquator, metersPerUnit, sphereRadius, useSLDEM);
            constexpr auto max = std::numeric_limits<uint8_t>::max();
            for(int i = 0; i < channelsPerPixel; ++i)
                outData[pixelPosInOutData + i] = std::clamp(normal[i]/2+0.5, 0.,1.) * max;
        }
    }
}

void cubeMapFindMaxAltitudes(const int16_t*const data, const ssize_t cubeMapSide, const QString& direction,
                             int16_t (&maxAltitudes)[MAX_ALT_TILES_PER_CUBE_SIDE][MAX_ALT_TILES_PER_CUBE_SIDE])
{
    std::cerr << "Finding maximum altitudes in sector '" << direction.toStdString() << "'...\n";
    if(cubeMapSide % MAX_ALT_TILES_PER_CUBE_SIDE)
        throw std::runtime_error("Cube map side isn't a multiple of max-altitude tiles count per side");

    const auto tileSide = cubeMapSide / MAX_ALT_TILES_PER_CUBE_SIDE;
    const auto numColsInTile = tileSide;
    const auto numRowsInTile = tileSide;

    for(int j = 0; j < MAX_ALT_TILES_PER_CUBE_SIDE; ++j)
    {
        for(int i = 0; i < MAX_ALT_TILES_PER_CUBE_SIDE; ++i)
        {
            const int16_t*const tile = data + tileSide*(cubeMapSide*j + i);
            int16_t maxDataValue = 0;
            for(int rowIdxInTile = 0; rowIdxInTile < numRowsInTile; ++rowIdxInTile)
            {
                const auto currRowStart = tile + rowIdxInTile*cubeMapSide;
                const auto currValue = *std::max_element(currRowStart, currRowStart + numColsInTile);
                if(currValue > maxDataValue) maxDataValue = currValue;
            }
            maxAltitudes[i][j] = maxDataValue;
            std::cerr << " - Tile [" << i << ", " << j << "]: " << maxDataValue << " m\n";
        }
    }
}

std::unique_ptr<QFile> makeCubeMapFaceFile(const QString& outDir, const ssize_t cubeMapSide, const QString& direction)
{
    const auto filePath = outDir + "/cube-face-" + direction + ".dat";
    std::cerr << "Creating " << filePath.toStdString() << "...\n";
    std::unique_ptr<QFile> file(new QFile(filePath));
    if(!file->open(QFile::ReadWrite))
    {
        throw std::runtime_error(QString("Failed to open %1 for writing: %2")
                                    .arg(filePath).arg(file->errorString()).toStdString());
    }

    const CubeMapFileHeader header{.formatVersion = CUBEMAP_FORMAT_VERSION, .cubeMapSide = uint32_t(cubeMapSide), .maxAltitudes={}};
    if(file->write(reinterpret_cast<const char*>(&header), sizeof header) != sizeof header)
    {
        throw std::runtime_error(QString("Failed to write header to %1: %2")
                                    .arg(filePath).arg(file->errorString()).toStdString());
    }

    return file;
}

void finalizeCubeMapFace(QFile& file, const ssize_t cubeMapSide, const QString& direction)
{
    CubeMapFileHeader header;
    file.seek(0);
    if(file.read(reinterpret_cast<char*>(&header), sizeof header) != sizeof header)
        throw std::runtime_error(QString("Failed to read back file header: %1").arg(file.errorString()).toStdString());

    const int16_t* data = nullptr;
    const qint64 sizeToRead = cubeMapSide * cubeMapSide * sizeof data[0];
    data = reinterpret_cast<const int16_t*>(file.map(sizeof header, sizeToRead));
    cubeMapFindMaxAltitudes(data, cubeMapSide, direction, header.maxAltitudes);
    file.seek(0);
    if(file.write(reinterpret_cast<const char*>(&header), sizeof header) != sizeof header)
    {
        throw std::runtime_error(QString("Failed to write %1 face final header: %2")
                                    .arg(direction).arg(file.errorString()).toStdString());
    }

    if(!file.flush())
    {
        throw std::runtime_error(QString("Failed to write %1 face: %2")
                                    .arg(direction).arg(file.errorString()).toStdString());
    }
    std::cerr << "Face creation completed\n";
}

void saveCubeMapChunk(QFile& file, const int16_t*const data, const ssize_t dataSize)
{
    const auto written = file.write(reinterpret_cast<const char*>(data), dataSize);
    if(written != dataSize)
    {
        throw std::runtime_error(QString("Failed to write data: wrote %1 bytes instead of %2")
                                    .arg(written).arg(dataSize).toStdString());
    }
}

void generateAltitudeCubeMap(const QString& outDir, const std::vector<Tile>& heightMapTiles, const TIFFTile& southLDEM,
                             const double resolutionAtEquator, const double metersPerUnit, const bool useSLDEM)
{
    // We include the edges of the cube in each face, so that on sampling we can linearly
    // interpolate the data inside a single face, without the need to split fetches.
#define PROCESS_FACE(X_EXPR, Y_EXPR, Z_EXPR, DIRECTION)                                                         \
    do {                                                                                                        \
        const auto file = makeCubeMapFaceFile(outDir, cubeMapSide, DIRECTION);                                  \
        std::vector<Type> dataChunk(linesPerChunk * cubeMapSide);                                               \
        ssize_t j;                                                                                              \
        for(j = 0; j < cubeMapSide; ++j)                                                                        \
        {                                                                                                       \
            const auto v = double(j) / (cubeMapSide-1);                                                         \
            for(ssize_t i = 0; i < cubeMapSide; ++i)                                                            \
            {                                                                                                   \
                const auto u = double(i) / (cubeMapSide-1);                                                     \
                const auto pixelPosInData = (i + j*cubeMapSide);                                                \
                const auto x = X_EXPR;                                                                          \
                const auto y = Y_EXPR;                                                                          \
                const auto z = Z_EXPR;                                                                          \
                                                                                                                \
                const auto longitude = atan2(y,x);                                                              \
                const auto latitude = std::asin(z / std::sqrt(x*x+y*y+z*z));                                    \
                                                                                                                \
                const auto samp = sample(heightMapTiles, southLDEM, resolutionAtEquator,                        \
                                         longitude, latitude, useSLDEM);                                        \
                const auto value = std::lround(metersPerUnit*samp);                                             \
                if(value < dataTypeMin || value > dataTypeMax)                                                  \
                {                                                                                               \
                    throw std::runtime_error(QString("Found a value out of data type range: %1")                \
                                                .arg(value).toStdString());                                     \
                }                                                                                               \
                dataChunk[pixelPosInData % dataChunk.size()] = value;                                           \
            }                                                                                                   \
            if((j+1) % linesPerChunk == 0)                                                                      \
            {                                                                                                   \
                saveCubeMapChunk(*file, dataChunk.data(), linesPerChunk * cubeMapSide * sizeof dataChunk[0]);   \
                std::cerr << j+1 << " lines done out of " << cubeMapSide << "\n";                               \
            }                                                                                                   \
        }                                                                                                       \
        saveCubeMapChunk(*file, dataChunk.data(), (j % linesPerChunk) * cubeMapSide * sizeof dataChunk[0]);     \
        finalizeCubeMapFace(*file, cubeMapSide, DIRECTION);                                                     \
    } while(false)

    const auto sourcePixelsPerCubeMapSide = resolutionAtEquator/4;
    // The scale coefficient takes into account the fact that pixel coordinates are now not, say,
    // longitude, but tan(longitude), and at the border of the cube map side it's tan(PI/4)=1, so,
    // to preserve the sampling rate at the densest part of the original data set, we need to divide
    // the resolution by PI/4.
    const ssize_t cubeMapSide0 = std::ceil(4/M_PI * sourcePixelsPerCubeMapSide);
    // To simplify work with subtiles, round cube map side to a multiple equal to the number of subtiles.
    constexpr int pieces = MAX_ALT_TILES_PER_CUBE_SIDE;
    const ssize_t cubeMapSide = (cubeMapSide0 + (pieces-1)) / pieces * pieces;

    using Type = int16_t;
    constexpr auto dataTypeMin = std::numeric_limits<Type>::min();
    constexpr auto dataTypeMax = std::numeric_limits<Type>::max();
    constexpr int linesPerChunk = 1000;
    PROCESS_FACE(  u*2-1 ,   v*2-1 ,   1  , "north");
    PROCESS_FACE(  u*2-1 , -(v*2-1),  -1  , "south");
    PROCESS_FACE(  u*2-1 ,    -1   , v*2-1, "west");
    PROCESS_FACE(    1   ,   u*2-1 , v*2-1, "lon0");
    PROCESS_FACE(-(u*2-1),     1   , v*2-1, "east");
    PROCESS_FACE(   -1   , -(u*2-1), v*2-1, "lon180");
}

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} ldemDir outDir";
    s << R"(
Options:
 -h, --help                 This help message
 --sldem DIR                Use SLDEM2015 from DIR to fill in the stripe of 60°S..60°N
 --cubemap                  Generate a cube map of heights instead of a normal map for
                             later use in generation of a horizons map
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

#define TIFFCHECK(call) if(call < 0) throw std::runtime_error(std::string(__FILE__)+":"+std::to_string(__LINE__)+": "+#call+" returned an error")

}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    QString ldemDir;
    QString sldemDir;
    QString ldemSouthPartPath;
    QString outDir;
    QString imgFormat = "png";
    QString surveyTitle;
    QString description;
    QString frame;
    QString creator;
    QString hips_copyright;
    QString obs_copyright;
    QString hipsStatus = "public mirror clonable";
    bool cubeMap = false;
    bool useSLDEM = false;
    bool useSeparateSouthLDEM = false;

    int totalPositionalArgumentsFound = 0;
    for(int n = 1; n < argc; ++n)
    {
        if(argv[n][0]!='-')
        {
            // Must be a positional argument
            switch(totalPositionalArgumentsFound)
            {
            case 0:
                ldemDir = argv[n];
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

        if(arg == "--sldem")
        {
            GO_TO_PARAM();
            sldemDir = argv[n];
            useSLDEM = true;
        }
        else if(arg == "--ldem-south")
        {
            GO_TO_PARAM();
            ldemSouthPartPath = argv[n];
            useSeparateSouthLDEM = true;
        }
        else if(arg == "--format")
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
        else if(arg == "--cubemap")
        {
            cubeMap = true;
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(ldemDir.isEmpty())
    {
        std::cerr << "LDEM directory not specified\n";
        return usage(argv[0], 1);
    }
    if(outDir.isEmpty())
    {
        std::cerr << "Output directory not specified\n";
        return usage(argv[0], 1);
    }

    if(sldemDir.isEmpty())
    {
        std::cerr << "SLDEM2015 directory not specified, using only LDEM\n";
    }

    QString finalExt;
    if(!cubeMap)
    {
        if(imgFormat == "jpeg")
            finalExt = "jpg";
        else if(imgFormat == "webp" || imgFormat == "bmp" || imgFormat == "png" || imgFormat == "tiff")
            finalExt = imgFormat;
        else
            throw std::runtime_error("Unexpected output image format: " + imgFormat.toStdString());
    }

    TIFFTile southLDEM;
    if(useSeparateSouthLDEM)
    {
        const auto path = ldemSouthPartPath.toStdString();
        std::cerr << "Reading \"" << path << "\"...\n";
        TIFFSetWarningHandler(0);
        const auto tiff = XTIFFOpen(path.c_str(), "r");
        if(!tiff) throw std::runtime_error("Failed to open the south LDEM file using libtiff");
        const auto gtiff = GTIFNew(tiff);
        if(!gtiff) throw std::runtime_error("Failed to initialize geotiff library");

        uint32_t width = 0, height = 0, depth = 0;
        tdir_t bestDir = 0;
        for(int currDir = 0;; ++currDir)
        {
            uint32_t currWidth, currHeight;
            uint16_t currDepth;
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &currWidth));
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &currHeight));
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &currDepth));
            std::cerr << "Found TIFF page with size " << currWidth << u8"×" << currHeight << "...\n";
            if(currWidth > width && currHeight >= height)
            {
                bestDir = TIFFCurrentDirectory(tiff);
                width = currWidth;
                height = currHeight;
                depth = currDepth;
            }
            if(!TIFFReadDirectory(tiff)) break;
        }
        std::cerr << "Using the largest page, size: " << width << u8"×" << height << ", depth: " << depth << "\n";
        TIFFSetDirectory(tiff, bestDir);
        // Make sure this is the directory we wanted
        {
            uint32_t currWidth, currHeight;
            uint16_t currDepth;
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &currWidth));
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &currHeight));
            TIFFCHECK(TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &currDepth));
            if(currWidth != width || currHeight != height || currDepth != depth)
                throw std::runtime_error("Failed to load the desired TIFF directory, loaded something strange");
        }
        uint16_t format = 0;
        TIFFCHECK(TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &format));
        constexpr uint16_t expectedFormat = SAMPLEFORMAT_IEEEFP;
        if(format != expectedFormat)
        {
            throw std::runtime_error("Unexpected TIFF data format: expected "+
                                     std::to_string(expectedFormat)+", got "+std::to_string(format));
        }
        uint16_t samplesPerPixel = 0;
        TIFFCHECK(TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel));
        if(samplesPerPixel != 1)
        {
            throw std::runtime_error("Unexpected number of samples per pixel in south LDEM TIFF "
                                     "file: expected 1, got "+std::to_string(samplesPerPixel));
        }
        if(!TIFFIsTiled(tiff))
            throw std::runtime_error("South LDEM TIFF file isn't tiled. This is unexpected and unsupported.");

        GTIFDefn defn;
        if(!GTIFGetDefn(gtiff, &defn))
            throw std::runtime_error("Failed to get geotiff definition");
        if(defn.Model != ModelTypeProjected)
            throw std::runtime_error("Unexpected geotiff model type (expected projected)");

        const auto projection = GTIFGetProj4Defn(&defn);
        if(!projection) throw std::runtime_error("Failed to get proj4 projection");
        southLDEM.ctx = proj_context_create();
        if(!southLDEM.ctx) throw std::runtime_error("Failed to create proj4 context");
        std::string lonLat = "+proj=longlat ";
        switch(defn.Ellipsoid)
        {
        case Ellipse_WGS_84:
            lonLat = "+ellps=WGS84";
            break;
        case Ellipse_Clarke_1866:
            lonLat = "+ellps=clrk66";
            break;
        case Ellipse_Clarke_1880:
            lonLat = "+ellps=clrk80";
            break;
        case Ellipse_GRS_1980:
            lonLat = "+ellps=GRS80";
            break;
        default:
            if(defn.SemiMajor && defn.SemiMinor)
            {
                lonLat += "+a=" + std::to_string(defn.SemiMajor) +
                          " +b=" + std::to_string(defn.SemiMinor);
            }
        }

        southLDEM.pj = proj_create_crs_to_crs(southLDEM.ctx, lonLat.c_str(), projection, nullptr);
        GTIFFreeMemory(projection);
        if(!southLDEM.pj) throw std::runtime_error("Failed to create CRS-to-CRS");

        southLDEM.width = width;
        southLDEM.height = height;
        southLDEM.gtiff = gtiff;
        southLDEM.defn = defn;
        southLDEM.data.reset(new float[width*height]);
        uint32_t tileWidth = 0, tileHeight = 0;
        TIFFCHECK(TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tileWidth));
        TIFFCHECK(TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tileHeight));
        std::cerr << "Tile size: " << tileWidth << u8"×" << tileHeight << "\n";
        std::cerr << "Decoding the image... ";
        std::unique_ptr<float[]> tile(new float[tileWidth*tileHeight]);
        for(uint32_t j = 0; j < height; j += tileHeight)
        {
            for(uint32_t i = 0; i < width; i += tileWidth)
            {
                TIFFReadTile(tiff, tile.get(), i, j, 0, 0);
                const auto realTileWidth = i + tileWidth > width ? width - i : tileWidth;
                const auto realTileHeight = j + tileHeight > height ? height - j : tileHeight;
                for(uint32_t tileRow = 0; tileRow < realTileHeight; ++tileRow)
                    std::copy(tile.get() + tileRow * realTileWidth,
                              tile.get() + (tileRow+1) * realTileWidth,
                              southLDEM.data.get() + (j+tileRow) * width + i);
            }
        }
        std::cerr << "done\n";
    }

    const QString ldemSectors[] =
    {
        "LDEM_512_45N_90N_180_270", "LDEM_512_45N_90N_270_360", "LDEM_512_45N_90N_000_090", "LDEM_512_45N_90N_090_180",
        "LDEM_512_00N_45N_180_270", "LDEM_512_00N_45N_270_360", "LDEM_512_00N_45N_000_090", "LDEM_512_00N_45N_090_180",
        "LDEM_512_45S_00S_180_270", "LDEM_512_45S_00S_270_360", "LDEM_512_45S_00S_000_090", "LDEM_512_45S_00S_090_180",
        "LDEM_512_90S_45S_180_270", "LDEM_512_90S_45S_270_360", "LDEM_512_90S_45S_000_090", "LDEM_512_90S_45S_090_180",
    };
    const QString sldemSectors[] =
    {
        "SLDEM2015_512_30N_60N_180_225", "SLDEM2015_512_30N_60N_225_270", "SLDEM2015_512_30N_60N_270_315", "SLDEM2015_512_30N_60N_315_360",
        "SLDEM2015_512_30N_60N_000_045", "SLDEM2015_512_30N_60N_045_090", "SLDEM2015_512_30N_60N_090_135", "SLDEM2015_512_30N_60N_135_180",
        "SLDEM2015_512_00N_30N_180_225", "SLDEM2015_512_00N_30N_225_270", "SLDEM2015_512_00N_30N_270_315", "SLDEM2015_512_00N_30N_315_360",
        "SLDEM2015_512_00N_30N_000_045", "SLDEM2015_512_00N_30N_045_090", "SLDEM2015_512_00N_30N_090_135", "SLDEM2015_512_00N_30N_135_180",
        "SLDEM2015_512_30S_00S_180_225", "SLDEM2015_512_30S_00S_225_270", "SLDEM2015_512_30S_00S_270_315", "SLDEM2015_512_30S_00S_315_360",
        "SLDEM2015_512_30S_00S_000_045", "SLDEM2015_512_30S_00S_045_090", "SLDEM2015_512_30S_00S_090_135", "SLDEM2015_512_30S_00S_135_180",
        "SLDEM2015_512_60S_30S_180_225", "SLDEM2015_512_60S_30S_225_270", "SLDEM2015_512_60S_30S_270_315", "SLDEM2015_512_60S_30S_315_360",
        "SLDEM2015_512_60S_30S_000_045", "SLDEM2015_512_60S_30S_045_090", "SLDEM2015_512_60S_30S_090_135", "SLDEM2015_512_60S_30S_135_180",
    };
    std::vector<Tile> heightMapTiles;
    for(const QString& sector : ldemSectors)
        heightMapTiles.push_back(readTile(ldemDir, sector));
    if(useSLDEM)
    {
        for(const QString& sector : sldemSectors)
            heightMapTiles.push_back(readTile(sldemDir, sector));
    }
    for(unsigned i = LDEM_FIRST + 1; i <= LDEM_LAST; ++i)
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
    if(useSLDEM)
    {
        for(unsigned i = SLDEM_FIRST + 1; i <= SLDEM_LAST; ++i)
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
    }

    for(unsigned i = 1; i < heightMapTiles.size(); ++i)
    {
        if(heightMapTiles[i].metersPerUnit != heightMapTiles[i-1].metersPerUnit)
        {
            throw std::runtime_error(QString(u8"Unit scaling factors of sectors %1 and %2 don't match: %3 km/unit vs %4 km/unit")
                                        .arg(heightMapTiles[i].sector).arg(heightMapTiles[i-1].sector)
                                        .arg(heightMapTiles[i].metersPerUnit).arg(heightMapTiles[i-1].metersPerUnit)
                                        .toStdString());
        }
    }
    const double metersPerUnit = heightMapTiles[0].metersPerUnit;
    for(unsigned i = 1; i < heightMapTiles.size(); ++i)
    {
        if(heightMapTiles[i].sphereRadius != heightMapTiles[i-1].sphereRadius)
        {
            throw std::runtime_error(QString(u8"Sphere radii of sectors %1 and %2 don't match: %3 vs %4")
                                        .arg(heightMapTiles[i].sector).arg(heightMapTiles[i-1].sector)
                                        .arg(heightMapTiles[i].sphereRadius).arg(heightMapTiles[i-1].sphereRadius)
                                        .toStdString());
        }
    }
    const double sphereRadius = heightMapTiles[0].sphereRadius;

    double resolutionAtEquator = 0; // pixels per 360°
    for(int i = LDEM_FIRST; i < LDEM_FIRST + LDEM_HORIZ_TILES_COUNT; ++i)
        resolutionAtEquator += heightMapTiles[i].width;

    if(!QDir().mkpath(outDir))
        throw std::runtime_error("Failed to create directory \""+outDir.toStdString()+'"');

    if(cubeMap)
    {
        generateAltitudeCubeMap(outDir, heightMapTiles, southLDEM, resolutionAtEquator, metersPerUnit, useSLDEM);
        return 0;
    }

    const int orderMax = std::ceil(std::log2(resolutionAtEquator / (4. * HIPS_TILE_SIZE * M_SQRT2)));

    hipsSaveProperties(outDir, orderMax, imgFormat, surveyTitle, "planet-normal",
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
                     &southLDEM = std::as_const(southLDEM),
                     &currPix,resolutionAtEquator,metersPerUnit,sphereRadius,useSLDEM,
                     &numThreadsReportedFirstProgress,&itemsDone,&itemsSkipped]()
        {
            constexpr int channelsPerPixel = 3;
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
                fillFace(orderMax, pix, heightMapTiles, southLDEM, channelsPerPixel, resolutionAtEquator,
                         metersPerUnit, sphereRadius, useSLDEM, data.data());
                const QImage out(data.data(), HIPS_TILE_SIZE, HIPS_TILE_SIZE, channelsPerPixel*HIPS_TILE_SIZE, QImage::Format_RGB888);
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
