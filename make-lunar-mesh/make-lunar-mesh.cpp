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
    static QImage in("/home/ruslan/Downloads/Moon/CGI Moon Kit displacement map - ldem_16_uint.tif");
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

void applyHeightMap(std::vector<ch_vertex>& vertices)
{
    for(auto& v : vertices)
    {
        const auto r0 = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
        const auto longitude = std::atan2(v.x, -v.y); // basically the same as atan(y,x)+90°, but in [-PI,PI]
        const auto latitude = std::acos(-v.z/r0) - 0.5*M_PI;
        const auto altitudeKm = moonHeightMapKmPerUnit*sample(moonHeightMapBits,
                                                              moonHeightMapWidth, moonHeightMapHeight,
                                                              moonHeightMapRowStride, longitude, latitude);
        const auto radius = (moonHeightMapBaseRadiusKm+altitudeKm)/AU;
        const auto scale = radius / r0;
        v.x *= scale;
        v.y *= scale;
        v.z *= scale;
    }
}

template<typename Rep, typename Per>
double toSeconds(const std::chrono::duration<Rep,Per> d)
{
    using namespace std::chrono;
    return duration_cast<microseconds>(d).count() * 1e-6;
}

std::vector<ch_vertex> combine(std::vector<ch_vertex> const& verticesCoarse,
                               std::vector<ch_vertex> const& verticesFine)
{
    std::vector<ch_vertex> vertices;
    vertices.reserve(verticesFine.size()); // This is an upper limit on the size

    const double fineYmax = std::sin(10*M_PI/180);
    const double fineYmin = -fineYmax;

    for(const auto& v : verticesCoarse)
        if(v.y > fineYmax)
            vertices.push_back(v);

    for(const auto& v : verticesFine)
        if(fineYmin <= v.y && v.y <= fineYmax)
            vertices.push_back(v);

    for(const auto& v : verticesCoarse)
        if(v.y < fineYmin)
            vertices.push_back(v);

    return vertices;
}

std::vector<ch_vertex> loadVertices(std::string const& path)
{
    std::ifstream file(path, std::ios_base::binary);
    if(!file) throw std::runtime_error("Failed to open input file \""+path+"\"");

    file.seekg(0, std::ios_base::end);
    const auto size = file.tellg();
    file.seekg(0);
    if(!size) throw std::runtime_error("Empty input file \""+path+"\"");

    static_assert(sizeof(ch_vertex)==3*sizeof(double));
    if(size%(sizeof(ch_vertex)))
        throw std::runtime_error("File size is not a multiple of 3×float64 for \""+path+"\"");

    std::vector<ch_vertex> vertices(size/sizeof(ch_vertex));
    file.read(reinterpret_cast<char*>(vertices.data()), size);
    if(!file) throw std::runtime_error("Failed to read input file \""+path+"\"");

    return vertices;
}

std::vector<int> loadIndices(std::string const& path)
{
    std::ifstream file(path, std::ios_base::binary);
    if(!file) throw std::runtime_error("Failed to open input file \""+path+"\"");

    file.seekg(0, std::ios_base::end);
    const auto size = file.tellg();
    file.seekg(0);
    if(!size) throw std::runtime_error("Empty input file \""+path+"\"");

    const auto numIndices = size/sizeof(uint64_t);
    if(numIndices % 3) throw std::runtime_error("Number of indices must be a multiple of 3");

    std::vector<uint64_t> data(numIndices);
    file.read(reinterpret_cast<char*>(data.data()), size);
    std::vector<int> indices(data.size());
    std::transform(data.begin(), data.end(), indices.begin(),
                   [](const uint64_t x){ return int(x); });

    return indices;
}

bool saveToBin(std::string const& path, std::vector<ch_vertex> const& vertices,
               int const*const indices, const uint32_t indexCount)
{
    std::ofstream out(path, std::ios_base::binary);
    const uint32_t vertexCount = vertices.size();
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
        std::transform(indices, indices+indexCount, indicesToWrite.begin(),
                       [](const int x){ return uint32_t(x); });
        out.write(reinterpret_cast<const char*>(indicesToWrite.data()),
                  indicesToWrite.size()*sizeof indicesToWrite[0]);
    }
    out.close();
    return !!out;
}

int usage(const char*const argv0, const int ret)
{
    (ret ? std::cerr : std::cout) << argv0
        << " {--fine-vertices|-f} input-file-fine-vertices.bin"
           " {{--coarse-vertices|-c} input-file-coarse-vertices.bin | {--fine-indices|-i} input-file-fine-indices.bin}"
           " [-o output-file.obj] [-b output-file.bin]\n";
    return ret;
}

int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    std::string inFileNameFineVertices;
    std::string inFileNameFineIndices;
    std::string inFileNameCoarseVertices;
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

        if(arg == "--fine-vertices" || arg == "-f")
        {
            GO_TO_PARAM();
            inFileNameFineVertices = argv[n];
        }
        else if(arg == "--coarse-vertices" || arg == "-c")
        {
            GO_TO_PARAM();
            inFileNameCoarseVertices = argv[n];
        }
        else if(arg == "--fine-indices" || arg == "-i")
        {
            GO_TO_PARAM();
            inFileNameFineIndices = argv[n];
        }
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

    if(inFileNameCoarseVertices.empty() == inFileNameFineIndices.empty())
    {
        std::cerr << "Either coarse vertices must be present, or fine indices\n";
        return 1;
    }

    const bool coarseMeshPresent = !inFileNameCoarseVertices.empty();

    loadHeightMap();

    using namespace std::chrono;

    const auto t0 = steady_clock::now();
    const auto verticesCoarse = coarseMeshPresent ? loadVertices(inFileNameCoarseVertices) : std::vector<ch_vertex>{};
    const auto t1 = steady_clock::now();
    if(coarseMeshPresent)
        std::cerr << verticesCoarse.size() << " vertices loaded for coarse mesh in " << toSeconds(t1-t0) << " s\n";
    const auto t2 = steady_clock::now();
    const auto verticesFine = loadVertices(inFileNameFineVertices);
    const auto t3 = steady_clock::now();
    std::cerr << verticesFine.size() << " vertices loaded for fine mesh in " << toSeconds(t3-t2) << " s\n";

    std::vector<int> indicesLoaded;
    int* indices = nullptr;
    int numFaces = 0;
    std::vector<ch_vertex> combined;
    if(coarseMeshPresent)
    {
        const auto t4 = steady_clock::now();
        combined = combine(verticesCoarse, verticesFine);
        const auto t5 = steady_clock::now();
        std::cerr << combined.size() << " vertices in the combined mesh in " << toSeconds(t5-t4) << " s\n";

        std::cerr << "Creating faces...\n";
        const auto t6 = steady_clock::now();
        convhull_3d_build(combined.data(), combined.size(), &indices, &numFaces);
        const auto t7 = steady_clock::now();
        std::cerr << numFaces << " faces created in " << toSeconds(t7-t6) << " s\n";
    }
    else
    {
        combined = verticesFine;

        const auto t0 = steady_clock::now();
        indicesLoaded = loadIndices(inFileNameFineIndices);
        numFaces = indicesLoaded.size() / 3;
        indices = indicesLoaded.data();
        const auto t1 = steady_clock::now();
        std::cerr << indicesLoaded.size() << " indices loaded in " << toSeconds(t1-t0) << " s\n";
    }

    std::cerr << "Applying heights...\n";
    const auto t8 = steady_clock::now();
    applyHeightMap(combined);
    const auto t9 = steady_clock::now();
    std::cerr << "Heights applied in " << toSeconds(t9-t8) << " s\n";

    std::cerr << "Saving results to \"" << outBinFileName << "\"...\n";
    const auto t10 = steady_clock::now();
    const bool ok = saveToBin(outBinFileName, combined, indices, 3*numFaces);
    const auto t11 = steady_clock::now();
    if(ok)
        std::cerr << "Saved in " << toSeconds(t11-t10) << " s\n";
    else
        std::cerr << "Failed to save to \"" << outBinFileName << "\"\n";

    std::cerr << "Saving results to \"" << outObjFileName << "\"...\n";
    const auto t12 = steady_clock::now();
    for(auto& v : combined)
    {
        // Saving an OBJ file via ConvHull stores vertices in fixed-point X.6
        // format, so don't try using astronomical units as the unit of length
        v.x *= AU;
        v.y *= AU;
        v.z *= AU;
    }
    convhull_3d_export_obj(combined.data(), combined.size(), indices, numFaces, false, outObjFileName.c_str());
    const auto t13 = steady_clock::now();
    std::cerr << "Saved in " << toSeconds(t13-t12) << " s\n";
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
