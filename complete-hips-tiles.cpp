#include <cmath>
#include <thread>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <QCryptographicHash>
#include <QPainter>
#include <QImage>
#include <QDir>
#include "hips.hpp"
#include "timing.hpp"
#include "healpix.hpp"

constexpr int TILE_SIZE = 512;
const QString initialExt = "bmp"; // Should be uncompressed for saving and reloading speed

int usage(const char* argv0, const int ret)
{
    auto& s = ret ? std::cerr : std::cout;
    s << "Usage: " << argv0 << "{options...} hipsDir";
    s << R"(
Options:
 -h, --help                 This help message
)";
    return ret;
}

int main(int argc, char** argv)
try
{
    QString hipsDir;

    int totalPositionalArgumentsFound = 0;
    for(int n = 1; n < argc; ++n)
    {
        if(argv[n][0]!='-')
        {
            // Must be a positional argument
            switch(totalPositionalArgumentsFound)
            {
            case 0:
                hipsDir = argv[n];
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
        {
            return usage(argv[0], 0);
        }
        else
        {
            std::cerr << "Unknown switch " << argv[n] << "\n";
            return usage(argv[0], 1);
        }
    }
    if(hipsDir.isEmpty())
    {
        std::cerr << "HiPS directory not specified\n";
        return usage(argv[0], 1);
    }

    QFile props(hipsDir + "/properties");
    if(!props.open(QFile::ReadOnly))
        throw std::runtime_error("Failed to open HiPS properties file");
    const auto text = props.readAll().split('\n');
    int orderMax = -1;
    for(const auto& line : text)
    {
        if(line.startsWith("#") || line.trimmed().isEmpty()) continue;
        const auto keyval = line.split('=');
        if(keyval.size() != 2)
            throw std::runtime_error(QString("Bad properties line: \"%1\"").arg(QString(line)).toStdString());
        bool ok = false;
        if(keyval[0].trimmed() == "hips_order")
        {
            orderMax = keyval[1].trimmed().toUInt(&ok);
            if(!ok)
                throw std::runtime_error(QString("Failed to parse HiPS order \"%1\"").arg(orderMax).toStdString());
        }
    }
    if(orderMax <= 0)
        throw std::runtime_error("Failed to find HiPS order");

    generateLowerOrderTiles(orderMax, hipsDir, true);

    std::cerr << "Generating Allsky previews...\n";
    for(int order = 0; order <= 3; ++order)
        generateAllsky(order, hipsDir, 128);
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
