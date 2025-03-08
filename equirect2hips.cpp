#include <cmath>
#include <thread>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <QCryptographicHash>
#include <QPainter>
#include <QImage>
#include <QDir>
#include "timing.hpp"
#include "healpix.hpp"

template<typename T>
double fetch(T const* data, const ssize_t width, const ssize_t height,
             const size_t rowStride, const size_t channelsPerPixel,
             const size_t subpixelIndex, const ssize_t requestedX, const ssize_t requestedY)
{
    const auto x = (requestedX+width) % width;
    const auto y = std::clamp(requestedY, ssize_t(0), height-1);

    const auto rawValue = data[x*channelsPerPixel + subpixelIndex + y*rowStride];

    if(std::is_integral_v<T> && std::is_unsigned_v<T>)
        return static_cast<double>(rawValue) / std::numeric_limits<T>::max();
    if(std::is_integral_v<T> && std::is_signed_v<T>)
        return static_cast<double>(rawValue) / -std::numeric_limits<T>::min();
    if(std::is_floating_point_v<T>)
        return rawValue;
}

template<typename T>
double sample(T const* data, const size_t width, const size_t height,
              const size_t rowStride, const size_t channelsPerPixel,
              const size_t subpixelIndex, const double longitude, const double latitude)
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

    const auto pTopLeft     = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY);
    const auto pTopRight    = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY);
    const auto pBottomLeft  = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX  , floorY+1);
    const auto pBottomRight = fetch(data, width, height, rowStride, channelsPerPixel, subpixelIndex, floorX+1, floorY+1);

    const auto fracX = x - floorX;
    const auto fracY = y - floorY;

    const auto sampleLeft  = pTopLeft  + (pBottomLeft -pTopLeft) *fracY;
    const auto sampleRight = pTopRight + (pBottomRight-pTopRight)*fracY;

    return sampleLeft + (sampleRight-sampleLeft)*fracX;
}

double normalizeLatitude(double x)
{
    if(x < -M_PI/2) return x + 2*M_PI;
    if(x > M_PI/2) return x - 2*M_PI;
    return x;
}

double normalizeLongitude(double x)
{
    if(x < -M_PI) return x + 2*M_PI;
    if(x > M_PI) return x - 2*M_PI;
    return x;
}

constexpr int TILE_SIZE = 512;
const QString initialExt = "bmp"; // Should be uncompressed for saving and reloading speed

void fillFace(const int order, const int pix, const QImage& in, double* outData)
{
    const unsigned nside = 1u << order;
    int ix, iy, face;
    healpix_nest2xyf(nside, pix, &ix, &iy, &face);

    const auto data = in.bits();
    const auto width  = in.width();
    const auto height = in.height();
    const auto strideInBytes = in.bytesPerLine();
    if(strideInBytes % sizeof data[0])
    {
        throw std::runtime_error("Row stride of " + std::to_string(strideInBytes) +
                                 " bytes is not a multiple of a pixel, this is not supported");
    }
    const auto bitDepth = in.depth();
    assert(bitDepth % (8 * sizeof data[0]) == 0);
    const int channelsPerPixel = bitDepth / (8 * sizeof data[0]);
    const auto rowStride = strideInBytes / sizeof data[0];

    for(int y = 0; y < TILE_SIZE; ++y)
    {
        for(int x = 0; x < TILE_SIZE; ++x)
        {
            double theta, phi;
            healpix_xyf2ang(nside * TILE_SIZE,
                            ix * TILE_SIZE + x, iy * TILE_SIZE + y,
                            face, &theta, &phi);
            const double latitude = normalizeLatitude(M_PI/2 - theta);
            const double longitude = normalizeLongitude(M_PI+phi);
            assert(-M_PI <= longitude && longitude <= M_PI);
            assert(-M_PI/2 <= latitude && latitude <= M_PI/2);

            // HiPS coordinates are swapped, because they were designed from the "look from
            // sphere center" perspective, while we are making a look from outside a planet.
            const int i = y, j = x;

            const int pixelPosInOutData = (i + j*TILE_SIZE)*channelsPerPixel;
            for(int subpixelN = 0; subpixelN < channelsPerPixel; ++subpixelN)
            {
                outData[pixelPosInOutData + subpixelN] =
                    sample(data, width, height, rowStride, channelsPerPixel,
                           subpixelN, longitude, latitude);
            }
        }
    }
}

void createLowerOrderTile(const int order, const int pix, const QString& outDir, QImage::Format imgFormat)
{
    const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(order).arg((pix / 10000) * 10000);
    if(!QDir().mkpath(outPath))
        throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
    const auto outFileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(initialExt);
    const auto inFileTemplate = QString("%1/Norder%2/Dir%3/Npix%5.%4").arg(outDir).arg(order + 1).arg((pix / 10000) * 10000).arg(initialExt);

    QImage outImg(TILE_SIZE, TILE_SIZE, imgFormat);
    {
        QPainter p(&outImg);
        for(int j = 0; j < 2; ++j)
        {
            for(int i = 0; i < 2; ++i)
            {
                const int innerPix = pix * 4 + i * 2 + j;
                const auto path = inFileTemplate.arg(innerPix);
                QImage img(path);
                if(img.isNull())
                    throw std::runtime_error("Failed to open \""+path.toStdString()+'"');
                img = img.scaled(TILE_SIZE / 2, TILE_SIZE / 2, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
                p.drawImage(QPoint(TILE_SIZE / 2 * i, TILE_SIZE / 2 * j), img);
            }
        }
    }
    outImg.save(outFileName);
}

QString md5sum(QString const& filePath)
{
    QFile file(filePath);
    if(!file.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open \"" << filePath.toStdString() << "\" for reading to compute MD5 sum\n";
        goto error;
    }
    {
        QCryptographicHash hash(QCryptographicHash::Md5);
        if(hash.addData(&file))
            return hash.result().toHex();
    }
    std::cerr << "Failed to compute MD5 sum of \"" << filePath.toStdString() << "\"\n";
error:
    std::cerr << "The property \"source_md5\" will be left empty\n";
    return "";
}

QString formatDate()
{
    return QDateTime::currentDateTime().toUTC().toString("yyyy-MM-ddTHH:mmZ");
}

void handleProgressReporting(const size_t totalItemCount,
                             const std::chrono::time_point<std::chrono::steady_clock> startTime,
                             std::chrono::time_point<std::chrono::steady_clock>& time0,
                             std::atomic_int& numThreadsReportedFirstProgress,
                             size_t& itemsDoneInThisThreadAfterLastUpdate,
                             std::atomic<unsigned>& itemsDone)
{
    ++itemsDoneInThisThreadAfterLastUpdate;
    using namespace std::chrono;
    const auto time1 = steady_clock::now();
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

        itemsDone += itemsDoneInThisThreadAfterLastUpdate;
        itemsDoneInThisThreadAfterLastUpdate = 0;
        const auto progress = std::round(itemsDone * 10000. / totalItemCount)/10000.;
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

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} equirect-image.ext outDir";
    s << R"(
Options:
 -h, --help                 This help message
 --format FORMAT            Set hips_tile_format to FORMAT (only one value is supported)
 --title TITLE              Set obs_title to TITLE
 --type TYPE                Set type (a non-standard property) to TYPE
 -d, --desc DESCRIPTION     Set obs_description to DESCRIPTION
 --frame FRAME              Set hips_frame to FRAME
 -c, --copyright COPYRIGHT  Set obs_copyright to COPYRIGHT
 -s, --status STATUS        Set hips_status to STATUS
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    setenv("QT_IMAGEIO_MAXALLOC","2048", false);

    QString inFileName;
    QString outDir;
    QString imgFormat = "png";
    QString surveyTitle;
    QString surveyType;
    QString description;
    QString frame;
    QString copyright;
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
                inFileName = argv[n];
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
        else if(arg == "--type")
        {
            GO_TO_PARAM();
            surveyType = argv[n];
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
        else if(arg == "-c" || arg == "--copyright")
        {
            GO_TO_PARAM();
            copyright = argv[n];
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
    if(inFileName.isEmpty())
    {
        std::cerr << "Input file not specified\n";
        return usage(argv[0], 1);
    }
    if(outDir.isEmpty())
    {
        std::cerr << "Output directory not specified\n";
        return usage(argv[0], 1);
    }

    std::cerr << "Loading input image... ";
    const QImage img(inFileName);
    if(img.isNull())
        throw std::runtime_error("Failed to open input image");
    else
        std::cerr << "done\n";

    QString finalExt;
    if(imgFormat == "jpeg")
        finalExt = "jpg";
    else if(imgFormat == "tiff")
        finalExt = "tif";
    else if(imgFormat == "webp" || imgFormat == "bmp" || imgFormat == "png")
        finalExt = imgFormat;
    else
        throw std::runtime_error("Unexpected output image format: " + imgFormat.toStdString());

    const auto in = img.convertToFormat(img.isGrayscale()     ? QImage::Format_Grayscale8 :
                                        img.hasAlphaChannel() ? QImage::Format_ARGB32     :
                                                                QImage::Format_RGB32);
    if(in.isNull())
        throw std::runtime_error("Failed to convert input image to the working format");

    const int orderMax = std::ceil(std::log2(in.width() / (4 * TILE_SIZE * M_SQRT2)));
    // First create the tiles of the deepest level
    std::cerr << "Creating tiles of order " << orderMax << "...\n";
    {
        const int absolutePixMax = 12 * (1 << (2 * orderMax));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,orderMax,in,outDir,startTime,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            std::vector<double> data(in.depth() * TILE_SIZE * TILE_SIZE);

            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                fillFace(orderMax, pix, in, data.data());

                using OutType = uchar;
                std::vector<OutType> outBits;
                outBits.reserve(data.size());
                if(std::is_integral_v<OutType>)
                {
                    for(auto v : data)
                        outBits.push_back(std::clamp(v,0.,1.)*std::numeric_limits<OutType>::max());
                }
                else
                {
                    for(auto v : data)
                        outBits.push_back(v*std::numeric_limits<OutType>::max());
                }
                const int channelsPerPixel = in.depth() / (8 * sizeof data[0]);
                const QImage out(outBits.data(), TILE_SIZE, TILE_SIZE, TILE_SIZE*channelsPerPixel, in.format());
                const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(orderMax).arg((pix / 10000) * 10000);
                if(!QDir().mkpath(outPath))
                    throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
                const auto fileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(initialExt);
                if(!out.save(fileName))
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

    // Now merge the deepest-level tiles to create the ones with lower detail level
    for(int order = orderMax - 1; order >= 0; --order)
    {
        std::cerr << "Creating tiles of order " << order << "...\n";
        const int absolutePixMax = 12 * (1 << (2 * order));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,order,in,outDir,startTime,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                createLowerOrderTile(order, pix, outDir, in.format());
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

    std::cerr << "Converting tiles to " << imgFormat.toStdString() << "...\n";
    struct OrderAndPix
    {
        int order;
        int pix;
        OrderAndPix(int order, int pix) : order(order), pix(pix) {}
    };
    std::vector<OrderAndPix> jobsToDo;
    for(int order = 0; order <= orderMax; ++order)
    {
        const int pixMax = 12 * (1 << (2 * order));
        for(int pix = 0; pix < pixMax; ++pix)
            jobsToDo.emplace_back(order, pix);
    }
    std::atomic_int numThreadsReportedFirstProgress{0};
    std::atomic<unsigned> itemsDone{0};
    const auto startTime = std::chrono::steady_clock::now();
    auto work = [&jobsToDo=std::as_const(jobsToDo),outDir,startTime,finalExt,
                 &numThreadsReportedFirstProgress,&itemsDone](const int jobMin, const int jobMax)
     {
         auto time0 = std::chrono::steady_clock::now();
         size_t itemsDoneInThisThreadAfterLastUpdate = 0;
         for(int n = jobMin; n < jobMax; ++n)
         {
             const int order = jobsToDo[n].order;
             const int pix = jobsToDo[n].pix;
             const auto pathTemplate = QString("%1/Norder%2/Dir%3/Npix%4.%5").arg(outDir).arg(order).arg((pix / 10000) * 10000).arg(pix);
             const auto inFileName = pathTemplate.arg(initialExt);
             const auto outFileName = pathTemplate.arg(finalExt);
             QImage img(inFileName);
             if(img.isNull())
             {
                 std::cerr << "Failed to open recently saved file " << inFileName.toStdString() << "\n";
                 break;
             }
             if(!img.save(outFileName))
             {
                 std::cerr << "Failed to save output file "  << outFileName.toStdString() << "\n";
                 break;
             }
             if(!QFile(inFileName).remove())
                 std::cerr << "Warning: failed to remove " << inFileName.toStdString() << "\n";

             handleProgressReporting(jobsToDo.size(), startTime, time0, numThreadsReportedFirstProgress,
                                     itemsDoneInThisThreadAfterLastUpdate, itemsDone);
         }
     };
    const auto time0 = std::chrono::steady_clock::now();
    const auto numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for(size_t n = 0; n < numThreads; ++n)
    {
        const size_t jobMin = jobsToDo.size() / numThreads * n;
        const size_t jobMax = n+1 < numThreads ? jobsToDo.size() / numThreads * (n+1) : jobsToDo.size();
        threads.emplace_back(work, jobMin, jobMax);
    }
    for(auto& thread : threads)
        thread.join();
    auto time1 = std::chrono::steady_clock::now();
    std::cerr << "100% done in " << formatDeltaTime(time0, time1) << "\n";

    const auto propsPath = outDir+"/properties";
    QFile propsFile(propsPath);
    if(!propsFile.open(QFile::WriteOnly))
        throw std::runtime_error("Failed to open \""+propsPath.toStdString()+"\" for writing");
    {
        QTextStream props(&propsFile);
        props << "hips_order            = " << orderMax << "\n";
        props << "hips_order_min        = 0\n";
        props << "hips_tile_width       = " << TILE_SIZE << "\n";
        props << "hips_tile_format      = " << imgFormat << "\n";
        props << "dataproduct_type      = image\n";
        props << "obs_title             = " << surveyTitle << "\n";
        props << "hips_release_date     = " << formatDate() << "\n";
        props << "source_md5            = " << md5sum(inFileName) << "\n";
        if(!surveyType.isEmpty())
            props << "type                  = " << surveyType << "\n";
        if(!description.isEmpty())
            props << "obs_description       = " << description << "\n";
        props << "hips_frame            = " << frame << "\n";
        if(!copyright.isEmpty())
            props << "obs_copyright         = " << copyright << "\n";
        props << "hips_status           = " << hipsStatus << "\n";
        props << "hips_version          = 1.4\n";
    }
    if(!propsFile.flush())
        throw std::runtime_error("Failed to write properties file");

    if(frame.isEmpty())
        std::cerr << "Warning: hips_frame is not specified\n";
    if(surveyTitle.isEmpty())
        std::cerr << "Warning: obs_title is not specified\n";
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
