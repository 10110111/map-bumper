#include <cmath>
#include <chrono>
#include <vector>
#include <fstream>
#include <iostream>
#include <QImage>
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"

constexpr double AU = 149597870.691; // km
const double moonHeightMapKmPerUnit = 0.5e-3;
const double moonHeightMapBaseRadiusKm = 1727.4;
const uint16_t* moonHeightMapBits;
int moonHeightMapWidth, moonHeightMapHeight, moonHeightMapRowStride;

void loadHeightMap()
{
    static QImage in("/home/ruslan/Downloads/Celestial bodies and space/Moon/generated-for-stellarium/23k/mixed-height-map.png");
    if(in.isNull())
        throw std::runtime_error("Failed to open input file");
    if(!in.isGrayscale())
        throw std::runtime_error("Input image is not grayscale. Can't interpret it as a height map.");

    in = in.convertToFormat(QImage::Format_Grayscale16);
    if(in.isNull())
        throw std::runtime_error("Failed to convert input file to Format_Grayscale16.");

    moonHeightMapBits = reinterpret_cast<uint16_t*>(in.bits());
    moonHeightMapWidth  = in.width();
    moonHeightMapHeight = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof moonHeightMapBits[0])
    {
        throw std::runtime_error("Row stride of "+std::to_string(strideInBytes)+" bytes "
                                 "is not a multiple of a pixel, this is not supported");
    }
    moonHeightMapRowStride = strideInBytes / sizeof moonHeightMapBits[0];
    std::cerr << "Loaded height map: " << moonHeightMapWidth << "×" << moonHeightMapHeight << "\n";
}

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
    assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

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

ch_vertex applyHeightMap(ch_vertex const& v)
{
    const auto r0 = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    const auto longitude = std::atan2(v.x, -v.y); // basically the same as atan(y,x)+90°, but in [-PI,PI]
    const auto latitude = std::acos(-v.z/r0) - 0.5*M_PI;
    const auto altitudeKm = moonHeightMapKmPerUnit*sample(moonHeightMapBits,
                                                          moonHeightMapWidth, moonHeightMapHeight,
                                                          moonHeightMapRowStride, longitude, latitude);
    const auto radius = (moonHeightMapBaseRadiusKm+altitudeKm)/AU;
    const auto scale = radius / r0;
    return {v.x * scale,
            v.y * scale,
            v.z * scale};
}

struct Mesh
{
    std::vector<ch_vertex> vertices;
    std::vector<int> indices;
};

enum class Direction
{
    X,
    MinusX,
    Y,
    MinusY,
    Z,
    MinusZ,
};

ch_vertex rotateZPlaneToDir(ch_vertex const& v, Direction dir)
{
    const auto x = v.x, y = v.y, z = v.z;
    switch(dir)
    {
    case Direction::X:
        return {z, y, -x};
    case Direction::MinusX:
        return {-z, -y, -x};
    case Direction::Y:
        return {-y, z, -x};
    case Direction::MinusY:
        return {y, -z, -x};
    case Direction::Z:
        return v;
    case Direction::MinusZ:
        return {-x, y, -z};
    }
    throw std::logic_error("bad direction");
}

// The whole sphere is split into 6 sides corresponding to a circumscribed cube
// onto which the vertices are projected from the center. Each side is split to
// an N×N grid of cells.
// This function creates the specified single cell from the upper (z=1) side of
// the cube.
Mesh createCell(const int numCellsPerCubeSide, const int cellIndex,
                const int numQuadsPerCellLength, const Direction dir)
{
    Mesh mesh;

    const double cellI = cellIndex % numCellsPerCubeSide;
    const double cellJ = cellIndex / numCellsPerCubeSide;
    const double cellX = 2. * cellI / numCellsPerCubeSide - 1;
    const double cellY = 2. * cellJ / numCellsPerCubeSide - 1;

    const double numQuadsPerSide = numQuadsPerCellLength * numCellsPerCubeSide;

#define RECT_GRID 1
#define ISOS_TRIANG_GRID 2
#define GRID_CHOICE ISOS_TRIANG_GRID
    const double z = 1;
    for(double j = 0; j <= numQuadsPerCellLength; ++j)
    {
        const double y = cellY + 2 * j / numQuadsPerSide;
        for(double i = 0; i <= numQuadsPerCellLength; ++i)
        {
#if GRID_CHOICE == RECT_GRID
            const double x = cellX + 2 * i / numQuadsPerSide;
#elif GRID_CHOICE == ISOS_TRIANG_GRID // TODO: straighten the zigzag borders
            const double shift = std::lround(j) % 2 ? 0.5 : 0;
            const double x = cellX + 2 * (i + shift) / numQuadsPerSide;
#endif
            // Apply equi-angular cubemap transformation
            const double eaX = std::tan(M_PI/4 * x);
            const double eaY = std::tan(M_PI/4 * y);

            const auto rotated = rotateZPlaneToDir({eaX,eaY,z}, dir);

            const auto v = applyHeightMap(rotated);
            mesh.vertices.push_back(v);
        }
    }

    const int lineSize = numQuadsPerCellLength + 1;
    for(int j = 0; j < numQuadsPerCellLength; ++j)
    {
        for(int i = 0; i < numQuadsPerCellLength; ++i)
        {
#if GRID_CHOICE == RECT_GRID
            mesh.indices.push_back( j    * lineSize + i  );
            mesh.indices.push_back( j    * lineSize + i+1);
            mesh.indices.push_back((j+1) * lineSize + i  );
            mesh.indices.push_back((j+1) * lineSize + i  );
            mesh.indices.push_back( j    * lineSize + i+1);
            mesh.indices.push_back((j+1) * lineSize + i+1);
#elif GRID_CHOICE == ISOS_TRIANG_GRID // TODO: straighten the zigzag borders
            if(j % 2 == 0)
            {
                mesh.indices.push_back( j    * lineSize + i  );
                mesh.indices.push_back( j    * lineSize + i+1);
                mesh.indices.push_back((j+1) * lineSize + i  );
                mesh.indices.push_back((j+1) * lineSize + i  );
                mesh.indices.push_back( j    * lineSize + i+1);
                mesh.indices.push_back((j+1) * lineSize + i+1);
            }
            else
            {
                mesh.indices.push_back( j    * lineSize + i  );
                mesh.indices.push_back((j+1) * lineSize + i+1);
                mesh.indices.push_back((j+1) * lineSize + i  );
                mesh.indices.push_back( j    * lineSize + i  );
                mesh.indices.push_back( j    * lineSize + i+1);
                mesh.indices.push_back((j+1) * lineSize + i+1);
            }
#endif
        }
    }

    return mesh;
}

template<typename Rep, typename Per>
double toSeconds(const std::chrono::duration<Rep,Per> d)
{
    using namespace std::chrono;
    return duration_cast<microseconds>(d).count() * 1e-6;
}

bool saveToBin(std::string const& path, std::vector<ch_vertex> const& vertices,
               std::vector<int> const& indices)
{
    std::ofstream out(path, std::ios_base::binary);
    const uint32_t vertexCount = vertices.size();
    const uint32_t indexCount = indices.size();
    out.write(reinterpret_cast<const char*>(&vertexCount), sizeof vertexCount);
    out.write(reinterpret_cast<const char*>(&indexCount), sizeof indexCount);
    {
        std::vector<float> verticesToWrite(3*vertices.size());
        std::transform(&vertices[0].v[0], &vertices[0].v[0]+verticesToWrite.size(),
                       verticesToWrite.begin(), [](const double x){ return float(x); });
        out.write(reinterpret_cast<const char*>(verticesToWrite.data()),
                  verticesToWrite.size()*sizeof verticesToWrite[0]);
    }
    {
        std::vector<uint32_t> indicesToWrite(indexCount);
        std::transform(indices.begin(), indices.end(), indicesToWrite.begin(),
                       [](const int x){ return uint32_t(x); });
        out.write(reinterpret_cast<const char*>(indicesToWrite.data()),
                  indicesToWrite.size()*sizeof indicesToWrite[0]);
    }
    out.close();
    return !!out;
}

int usage(const char*const argv0, const int ret)
{
    (ret ? std::cerr : std::cout) << argv0 << " [-o output-file.obj] [-b output-file.bin]\n";
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    std::string outObjFileName;
    std::string outBinFileName;
    for(int n = 1; n < argc; ++n)
    {
        const std::string arg(argv[n]);
        if(arg == "--help" || arg == "-h")
            return usage(argv[0], 0);

#define GO_TO_PARAM()                                                       \
            ++n;                                                            \
            if(n == argc)                                                   \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return 1;                                                   \
            }                                                               \
            do{}while(0)

        else if(arg == "-o")
        {
            GO_TO_PARAM();
            outObjFileName = argv[n];
        }
        else if(arg == "-b")
        {
            GO_TO_PARAM();
            outBinFileName = argv[n];
        }
        else
        {
            return usage(argv[0], 1);
        }
    }

    using namespace std::chrono;

    std::cerr << "Loading height map...\n";
    const auto t0 = steady_clock::now();
    loadHeightMap();
    const auto t1 = steady_clock::now();
    std::cerr << "Loaded in " << toSeconds(t1-t0) << " s\n";

    auto mesh = createCell(8, 21, 1000, Direction::X);

    if(!outBinFileName.empty())
    {
        std::cerr << "Saving results to \"" << outBinFileName << "\"...\n";
        const auto t10 = steady_clock::now();
        const bool ok = saveToBin(outBinFileName, mesh.vertices, mesh.indices);
        const auto t11 = steady_clock::now();
        if(ok)
            std::cerr << "Saved in " << toSeconds(t11-t10) << " s\n";
        else
            std::cerr << "Failed to save to \"" << outBinFileName << "\"\n";
    }

    if(!outObjFileName.empty())
    {
        std::cerr << "Saving results to \"" << outObjFileName << "\"...\n";
        const auto t12 = steady_clock::now();
        for(auto& v : mesh.vertices)
        {
            // Saving an OBJ file via ConvHull stores vertices in fixed-point X.6
            // format, so don't try using astronomical units as the unit of length
            v.x *= AU;
            v.y *= AU;
            v.z *= AU;
        }
        convhull_3d_export_obj(mesh.vertices.data(), mesh.vertices.size(),
                               mesh.indices.data(), mesh.indices.size() / 3,
                               false, outObjFileName.c_str());
        const auto t13 = steady_clock::now();
        std::cerr << "Saved in " << toSeconds(t13-t12) << " s\n";
    }
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
