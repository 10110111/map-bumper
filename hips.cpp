#include "hips.hpp"
#include <thread>
#include <iostream>
#include <QDir>
#include <QFile>
#include <QPainter>
#include <QDateTime>
#include <QCryptographicHash>
#include "timing.hpp"

QString formatDate()
{
    return QDateTime::currentDateTime().toUTC().toString("yyyy-MM-ddTHH:mmZ");
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

void createLowerOrderTile(const int order, const int pix, const QString& outDir, const bool savePartialTiles)
{
    try
    {
        const auto outPath = QString("%1/Norder%2/Dir%3").arg(outDir).arg(order).arg((pix / 10000) * 10000);
        if(!QDir().mkpath(outPath))
            throw std::runtime_error("Failed to create directory \""+outPath.toStdString()+'"');
        const auto outFileName = QString("%1/Npix%4.%5").arg(outPath).arg(pix).arg(hipsInitialExt);
        const auto inFileTemplate = QString("%1/Norder%2/Dir%4/Npix%5.%3").arg(outDir).arg(order + 1).arg(hipsInitialExt);

        QImage outImg(HIPS_TILE_SIZE, HIPS_TILE_SIZE, QImage(inFileTemplate.arg(0).arg(0)).format());
        outImg.fill(QColor(0,0,0,0));
        {
            QPainter p(&outImg);
            for(int j = 0; j < 2; ++j)
            {
                for(int i = 0; i < 2; ++i)
                {
                    const int innerPix = pix * 4 + i * 2 + j;
                    const auto path = inFileTemplate.arg((innerPix / 10000) * 10000).arg(innerPix);
                    QImage img(path);
                    if(img.isNull())
                    {
                        if(savePartialTiles)
                        {
                            std::cerr << "Failed to open \""+path.toStdString()+"\"\n";
                            continue;
                        }
                        else
                        {
                            throw std::runtime_error("Failed to open \""+path.toStdString()+'"');
                        }
                    }
                    img = img.scaled(HIPS_TILE_SIZE / 2, HIPS_TILE_SIZE / 2, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
                    p.drawImage(QPoint(HIPS_TILE_SIZE / 2 * i, HIPS_TILE_SIZE / 2 * j), img);
                }
            }
        }
        outImg.save(outFileName, nullptr, 100);
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << "\n";
    }
}

void generateLowerOrderTiles(const int orderMax, const QString& outDir, const bool savePartialTiles)
{
    for(int order = orderMax - 1; order >= 0; --order)
    {
        std::cerr << "Creating tiles of order " << order << "...\n";
        const int absolutePixMax = 12 * (1 << (2 * order));
        std::atomic_int numThreadsReportedFirstProgress{0};
        std::atomic<unsigned> itemsDone{0};
        const auto startTime = std::chrono::steady_clock::now();
        auto work = [absolutePixMax,order,outDir,startTime,savePartialTiles,
                     &numThreadsReportedFirstProgress,&itemsDone](const int pixMin, const int pixMax)
        {
            auto time0 = std::chrono::steady_clock::now();
            size_t itemsDoneInThisThreadAfterLastUpdate = 0;
            for(int pix = pixMin; pix < pixMax; ++pix)
            {
                createLowerOrderTile(order, pix, outDir, savePartialTiles);
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

}

void generateAllsky(const int order, QString const& hipsDir, const int tileSize)
{
    const int pixMax = 12 * (1 << (2 * order));
    const int numTilesInRow = int(std::sqrt(pixMax));
    const int numTilesRows = (pixMax + numTilesInRow - 1) / numTilesInRow;
    const auto inFileTemplate = QString("%1/Norder%2/Dir%4/Npix%5.%3").arg(hipsDir).arg(order).arg(hipsInitialExt);
    QImage outImg(tileSize * numTilesInRow, tileSize * numTilesRows, QImage(inFileTemplate.arg(0).arg(0)).format());
    outImg.fill(QColor(0,0,0,0));
    {
        QPainter p(&outImg);
        for(int pix = 0; pix < pixMax; ++pix)
        {
            const auto path = inFileTemplate.arg((pix / 10000) * 10000).arg(pix);
            QImage tile(path);
            if(tile.isNull())
                throw std::runtime_error("Failed to open \""+path.toStdString()+'"');
            tile = tile.scaled(tileSize, tileSize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
            p.drawImage(QPoint(tileSize * (pix % numTilesInRow), tileSize * (pix / numTilesInRow)), tile);
        }
    }
    const auto outFileName = QString("%1/Norder%2/Allsky.%3").arg(hipsDir).arg(order).arg(hipsInitialExt);
    outImg.save(outFileName, nullptr, 100);
}

void convertTiles(const QString& finalExt, const QString& formatName, const int orderMax, const QString& outDir)
{
    if(finalExt == hipsInitialExt) return;

    std::cerr << "Converting tiles to " << formatName.toStdString() << "...\n";
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
             const auto inFileName = pathTemplate.arg(hipsInitialExt);
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
}

void hipsSaveProperties(QString const& outDir, const int orderMax, QString const& imgFormat, QString const& surveyTitle,
                        QString const& surveyType, QString const& description, QString const& frame, QString const& obs_copyright,
                        QString const& hips_copyright, QString const& creator, QString const& hipsStatus)
{
    const auto propsPath = outDir+"/properties";
    QFile propsFile(propsPath);
    if(!propsFile.open(QFile::WriteOnly))
        throw std::runtime_error("Failed to open \""+propsPath.toStdString()+"\" for writing");
    {
        QTextStream props(&propsFile);
        props << "hips_order            = " << orderMax << "\n";
        props << "hips_order_min        = 0\n";
        props << "hips_tile_width       = " << HIPS_TILE_SIZE << "\n";
        props << "hips_tile_format      = " << imgFormat << "\n";
        props << "dataproduct_type      = image\n";
        props << "obs_title             = " << surveyTitle << "\n";
        props << "hips_release_date     = " << formatDate() << "\n";
        if(!surveyType.isEmpty())
            props << "type                  = " << surveyType << "\n";
        if(!description.isEmpty())
            props << "obs_description       = " << description << "\n";
        props << "hips_frame            = " << frame << "\n";
        if(!obs_copyright.isEmpty())
            props << "obs_copyright         = " << obs_copyright << "\n";
        if(!hips_copyright.isEmpty())
            props << "hips_copyright        = " << hips_copyright << "\n";
        if(!creator.isEmpty())
            props << "hips_creator          = " << creator << "\n";
        props << "hips_version          = 1.4\n";
        props << "hips_status           = " << hipsStatus << "\n";
    }
    if(!propsFile.flush())
        throw std::runtime_error("Failed to write properties file");

    if(frame.isEmpty())
        std::cerr << "Warning: hips_frame is not specified\n";
    if(surveyTitle.isEmpty())
        std::cerr << "Warning: obs_title is not specified\n";
}
