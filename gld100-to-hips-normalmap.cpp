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

struct Tile
{
    QString sector;
    double sphereRadius;
    double metersPerUnit;
    double lineProjectionOffset;
    double sampleProjectionOffset;
    double mapScale;
    std::unique_ptr<QFile> file;
    const int16_t* data;
    ssize_t width, height;
};

enum SectorOffset
{
    P900N_256P,
    P900N,
    N2250, N3150, N0450, N1350,
    S2250, S3150, S0450, S1350,
    P900S,
    P900S_256P,
};

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

std::pair<double/*i*/,double/*j*/> latLonToStereoPoint(const double longitude, const double latitude, const double sphereRadius,
                                                       const double mapScale, const double lineProjectionOffset,
                                                       const double sampleProjectionOffset, bool north)
{
    // Following the description in DSMAP.CAT
    using namespace std;
    double x, y;
    if(north)
    {
        const auto tanL = tan(M_PI/4 - latitude/2);
        x =  2*sphereRadius*sin(longitude)*tanL;
        y =  2*sphereRadius*cos(longitude)*tanL;
    }
    else
    {
        const auto tanL = tan(M_PI/4 + latitude/2);
        x =  2*sphereRadius*sin(longitude)*tanL;
        y = -2*sphereRadius*cos(longitude)*tanL;
    }
    const auto i = sampleProjectionOffset + x / mapScale;
    const auto j = lineProjectionOffset + y / mapScale;
    return {i,j};
}

double fetch(std::vector<Tile> const& data, const ssize_t requestedX, const ssize_t requestedY, const ssize_t resolutionAtEquator)
{
    const auto totalWidth = resolutionAtEquator;
    const auto totalHeight = totalWidth / 2;
    const auto x = (requestedX+totalWidth) % totalWidth;
    auto y = std::clamp(requestedY, ssize_t(0), totalHeight-1);

    const auto polarSectorHeight = data[N2250].height / 2;
    const auto northSectorBottomY = polarSectorHeight;
    const auto southSectorTopY = northSectorBottomY + 2 * data[N2250].height;
    const auto equatSectorWidth = data[N2250].width;
    const auto equatSectorHeight = data[N2250].height;

    if(y < northSectorBottomY)
    {
        std::cerr << "fetch: warning: y got too far north: " << y << "\n";
        y = northSectorBottomY;
    }
    else if(y >= southSectorTopY)
    {
        std::cerr << "fetch: warning: y got too far south: " << y << "\n";
        y = southSectorTopY;
    }

    if(x < totalWidth / 2)
    {
        if(x < data[N2250].width)
        {
            // E300.2250
            if(y < totalHeight / 2)
            {
                // E300N2250
                return data[N2250].data[(y-northSectorBottomY)*equatSectorWidth + x];
            }
            else
            {
                // E300S2250
                return data[S2250].data[(y-(northSectorBottomY+equatSectorHeight))*equatSectorWidth + x];
            }
        }
        else
        {
            // E300.3150
            if(y < totalHeight / 2)
            {
                // E300N3150
                return data[N3150].data[(y-northSectorBottomY)*equatSectorWidth + (x - equatSectorWidth)];
            }
            else
            {
                // E300S3150
                return data[S3150].data[(y-(northSectorBottomY+equatSectorHeight))*equatSectorWidth + (x - equatSectorWidth)];
            }
        }
    }
    else
    {
        if(x < 3 * data[N2250].width)
        {
            // E300.0450
            if(y < totalHeight / 2)
            {
                // E300N0450
                return data[N0450].data[(y-northSectorBottomY)*equatSectorWidth + (x - 2*equatSectorWidth)];
            }
            else
            {
                // E300S0450
                return data[S0450].data[(y-(northSectorBottomY+equatSectorHeight))*equatSectorWidth + (x - 2*equatSectorWidth)];
            }
        }
        else
        {
            // E300.1350
            if(y < totalHeight / 2)
            {
                // E300N1350
                return data[N1350].data[(y-northSectorBottomY)*equatSectorWidth + (x - 3*equatSectorWidth)];
            }
            else
            {
                // E300S1350
                return data[S1350].data[(y-(northSectorBottomY+equatSectorHeight))*equatSectorWidth + (x - 3*equatSectorWidth)];
            }
        }
    }
}

double fetchFromRect(int16_t const* data, const ssize_t width, const ssize_t height,
                     const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = std::clamp(requestedX, ssize_t(0), width-1);
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);
    return data[y*width + x];
}

bool isBad(const int16_t value)
{
    return value == -32768;
}

double samplePolarSector(int16_t const* data, const size_t width, const size_t height,
                         const double x, const double y)
{
    const auto floorX = std::floor(x);
    const auto floorY = std::floor(y);

    const auto pTopLeft     = fetchFromRect(data, width, height, floorX  , floorY);
    const auto pTopRight    = fetchFromRect(data, width, height, floorX+1, floorY);
    const auto pBottomLeft  = fetchFromRect(data, width, height, floorX  , floorY+1);
    const auto pBottomRight = fetchFromRect(data, width, height, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = isBad(pTopLeft) ? pBottomLeft : isBad(pBottomLeft) ? pTopLeft : pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = isBad(pTopRight)? pBottomRight: isBad(pBottomRight)? pTopRight: pTopRight + (pBottomRight-pTopRight)*fracY;
    if(isBad(sampleLeft)) return sampleRight;
    if(isBad(sampleRight)) return sampleLeft;


    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

double sampleNorthSector(std::vector<Tile> const& data, const double longitude, const double latitude)
{
    const auto& sec = data[P900N];
    const auto [i,j] = latLonToStereoPoint(longitude, latitude, sec.sphereRadius, sec.mapScale,
                                           sec.lineProjectionOffset, sec.sampleProjectionOffset, true);
    const auto value = samplePolarSector(sec.data, sec.width, sec.height, i, j);

    if(isBad(value))
    {
        // Fill the hole at the pole
        const auto& sec = data[P900N_256P];
        const auto [i,j] = latLonToStereoPoint(longitude, latitude, sec.sphereRadius, sec.mapScale,
                                               sec.lineProjectionOffset, sec.sampleProjectionOffset, true);
        return samplePolarSector(sec.data, sec.width, sec.height, i, j);
    }

    return value;
}

double sampleSouthSector(std::vector<Tile> const& data, const double longitude, const double latitude)
{
    const auto& sec = data[P900S];
    const auto [i,j] = latLonToStereoPoint(longitude, latitude, sec.sphereRadius, sec.mapScale,
                                           sec.lineProjectionOffset, sec.sampleProjectionOffset, false);
    const auto value = samplePolarSector(sec.data, sec.width, sec.height, i, j);

    if(isBad(value))
    {
        // Fill the hole at the pole
        const auto& sec = data[P900S_256P];
        const auto [i,j] = latLonToStereoPoint(longitude, latitude, sec.sphereRadius, sec.mapScale,
                                               sec.lineProjectionOffset, sec.sampleProjectionOffset, false);
        return samplePolarSector(sec.data, sec.width, sec.height, i, j);
    }

    return value;
}

double sample(std::vector<Tile> const& data, const double resolutionAtEquator, double longitude, double latitude)
{
    if(latitude > M_PI/2 || latitude < -M_PI/2)
    {
        latitude = -latitude;
        longitude += M_PI;
    }

    const auto deltaLon = 2.*M_PI/resolutionAtEquator;
    const auto firstLon = (1.-resolutionAtEquator)/2. * deltaLon;

    const auto x = (longitude - firstLon) / deltaLon;
    const auto floorX = std::floor(x);

    const auto deltaLat = -M_PI/(resolutionAtEquator/2-1);
    const auto firstLat = (1.-resolutionAtEquator/2)/2. * deltaLat;

    const auto y = (latitude - firstLat) / deltaLat;
    const auto floorY = std::floor(y);

    const auto polarSectorHeight = data[N2250].height / 2;
    const auto northSectorBottomY = polarSectorHeight;
    if(floorY < northSectorBottomY)
        return sampleNorthSector(data, longitude, latitude);
    const auto southSectorTopY = northSectorBottomY + 2 * data[N2250].height;
    if(floorY+1 >= southSectorTopY)
        return sampleSouthSector(data, longitude, latitude);

    const auto pTopLeft     = fetch(data, floorX  , floorY  , resolutionAtEquator);
    const auto pTopRight    = fetch(data, floorX+1, floorY  , resolutionAtEquator);
    const auto pBottomLeft  = fetch(data, floorX  , floorY+1, resolutionAtEquator);
    const auto pBottomRight = fetch(data, floorX+1, floorY+1, resolutionAtEquator);

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
    const auto filename = inDir + "/WAC_GLD100_" + sector + ".IMG";
    tile.file.reset(new QFile(filename));
    if(!tile.file->open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << filename.toStdString()
                  << " for reading: "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    QByteArray header(1024, '\0');
    if(tile.file->read(header.data(), header.size()) != header.size())
    {
        std::cerr << "Failed to read header of "
                  << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    const auto headerSize = QRegularExpression("RECORD_BYTES\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    if(headerSize == 0)
    {
        std::cerr << "Failed to determine header size\n";
        return {};
    }
    header = QByteArray(headerSize, '\0');
    tile.file->seek(0);
    if(tile.file->read(header.data(), header.size()) != header.size())
    {
        std::cerr << "Failed to read header of "
                  << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }

    tile.width = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    tile.height= QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    if(tile.width==0 || tile.height==0)
    {
        std::cerr << "Failed to find image dimensions\n";
        return {};
    }
    std::cerr << "Image dimensions for sector " << sector.toStdString()
              << ": " << tile.width << u8"×" << tile.height << "\n";
    const auto radiusA = QRegularExpression(R"(\n\s*A_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<KM>\s*\n)").match(header).captured(1).toDouble();
    const auto radiusB = QRegularExpression(R"(\n\s*B_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<KM>\s*\n)").match(header).captured(1).toDouble();
    const auto radiusC = QRegularExpression(R"(\n\s*C_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<KM>\s*\n)").match(header).captured(1).toDouble();
    if(radiusA <= 0 || radiusA != radiusB || radiusB != radiusC)
        throw std::runtime_error(QString("Bad reference sphere radii found. Expected three equal positive radii, but found %1, %2, %3")
                                    .arg(radiusA).arg(radiusB).arg(radiusC).toStdString());
    tile.sphereRadius = 1000 * radiusA;

    const auto sampleOffset = QRegularExpression("\\bOFFSET\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toDouble();
    const auto sampleScalingFactor = QRegularExpression("\\bSCALING_FACTOR\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toDouble();
    if(sampleOffset != 0)
        throw std::runtime_error(QString("Unexpected sample offset: %1").arg(sampleOffset).toStdString());
    tile.metersPerUnit = sampleScalingFactor;
    const auto sampleBits = QRegularExpression("\\bSAMPLE_BITS\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
    if(sampleBits != 16)
        throw std::runtime_error("Unexpected sample size of "+std::to_string(sampleBits)+" bits");

    tile.lineProjectionOffset = QRegularExpression(R"(\n\s*LINE_PROJECTION_OFFSET\s*=\s*([0-9.]+)\s*<PIXEL>\s*\n)").match(header).captured(1).toDouble();
    tile.sampleProjectionOffset = QRegularExpression(R"(\n\s*SAMPLE_PROJECTION_OFFSET\s*=\s*([0-9.]+)\s*<PIXEL>\s*\n)").match(header).captured(1).toDouble();
    tile.mapScale = QRegularExpression(R"(\n\s*MAP_SCALE\s*=\s*([0-9.]+)\s*<METERS/PIXEL>\s*\n)").match(header).captured(1).toDouble();
    if(tile.sector.startsWith('P'))
    {
        if(!tile.lineProjectionOffset || !tile.sampleProjectionOffset || !tile.mapScale)
            throw std::runtime_error(QString("Bad projection parameters: line projection offset = %1, sample projection offset = %2, map scale = %3")
                                        .arg(tile.lineProjectionOffset).arg(tile.sampleProjectionOffset).arg(tile.mapScale).toStdString());
    }

    const qint64 sizeToRead = tile.width * tile.height * sizeof tile.data[0];
    tile.data = reinterpret_cast<const int16_t*>(tile.file->map(header.size(), sizeToRead));
    if(!tile.data)
    {
        std::cerr << "Failed to read " << filename.toStdString() << ": "
                  << tile.file->errorString().toStdString() << "\n";
        return {};
    }
    return tile;
}

glm::dvec3 computeNormal(const double centerLon, const double centerLat,
                         const std::vector<Tile>& heightMapTiles,
                         const double resolutionAtEquator,
                         const double metersPerUnit, const double sphereRadius)
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

    const double eastHeight = metersPerUnit*sample(heightMapTiles, resolutionAtEquator, eastLon, eastLat);
    const double eastRadius = sphereRadius + eastHeight;

    const double westHeight = metersPerUnit*sample(heightMapTiles, resolutionAtEquator, westLon, westLat);
    const double westRadius = sphereRadius + westHeight;

    const double northHeight = metersPerUnit*sample(heightMapTiles, resolutionAtEquator, northLon, northLat);
    const double northRadius = sphereRadius + northHeight;

    const double southHeight = metersPerUnit*sample(heightMapTiles, resolutionAtEquator, southLon, southLat);
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
              const int channelsPerPixel, const double resolutionAtEquator,
              const double metersPerUnit, const double sphereRadius, uint8_t* outData)
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
            const auto normal = computeNormal(longitude, latitude, heightMapTiles,
                                              resolutionAtEquator, metersPerUnit, sphereRadius);
            constexpr auto max = std::numeric_limits<uint8_t>::max();
            for(int i = 0; i < channelsPerPixel; ++i)
                outData[pixelPosInOutData + i] = std::clamp(normal[i]/2+0.5, 0.,1.) * max;
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

    const QString sectors[] = {
                               "P900N0000_256P",
                               "P900N0000_100M",
                               "E300N2250_100M", "E300N3150_100M", "E300N0450_100M", "E300N1350_100M",
                               "E300S2250_100M", "E300S3150_100M", "E300S0450_100M", "E300S1350_100M",
                               "P900S0000_100M",
                               "P900S0000_256P",
                              };
    std::vector<Tile> heightMapTiles;
    for(const QString& sector : sectors)
        heightMapTiles.push_back(readTile(inDir, sector));
    for(unsigned i = N2250+1; i <= S1350; ++i)
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
    for(unsigned i = N2250; i <= N1350; ++i)
        resolutionAtEquator += heightMapTiles[i].width;

    const int orderMax = std::ceil(std::log2(resolutionAtEquator / (4. * HIPS_TILE_SIZE * M_SQRT2)));

    // First create the tiles of the deepest level
    std::cerr << "Creating tiles of order " << orderMax << "...\n";
    {
        const int absolutePixMax = 12 * (1 << (2 * orderMax));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,orderMax,outDir,startTime,&heightMapTiles = std::as_const(heightMapTiles),
                     resolutionAtEquator,metersPerUnit,sphereRadius,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            constexpr int channelsPerPixel = 3;
            std::vector<uint8_t> data(HIPS_TILE_SIZE * HIPS_TILE_SIZE * channelsPerPixel);

            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                fillFace(orderMax, pix, heightMapTiles, channelsPerPixel, resolutionAtEquator, metersPerUnit, sphereRadius, data.data());
                const QImage out(data.data(), HIPS_TILE_SIZE, HIPS_TILE_SIZE, channelsPerPixel*HIPS_TILE_SIZE, QImage::Format_RGB888);
                const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(orderMax).arg((pix / 10000) * 10000);
                if(!QDir().mkpath(outPath))
                    throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
                const auto fileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(hipsInitialExt);
                if(!out.save(fileName, nullptr, 100))
                    throw std::runtime_error("Failed to save output file " + fileName.toStdString());

                handleProgressReporting(absolutePixMax, startTime, time0, numThreadsReportedFirstProgress,
                                        itemsDoneInThisThreadAfterLastUpdate, itemsDone);
            }
        };
        const auto time0 = std::chrono::steady_clock::now();
        const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        for(size_t n = 0; n < numThreads; ++n)
        {
            const size_t pixMin = absolutePixMax / numThreads * n;
            const size_t pixMax = n+1 < numThreads ? absolutePixMax / numThreads * (n+1) : absolutePixMax;
            threads.emplace_back(work, pixMin, pixMax);
        }
        for(auto& thread : threads)
            thread.join();
        auto time1 = std::chrono::steady_clock::now();
        std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";
    }

    generateLowerOrderTiles(orderMax, outDir);
    convertTiles(finalExt, imgFormat, orderMax, outDir);
    hipsSaveProperties(outDir, orderMax, imgFormat, surveyTitle, "planet-normal", description,
                       frame, obs_copyright, hips_copyright, creator, hipsStatus);
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}

