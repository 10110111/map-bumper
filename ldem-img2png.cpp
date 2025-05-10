#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <QRegularExpression>
#include <QString>
#include <QImage>
#include <QFile>

constexpr double DEGREE = M_PI/180;

int usage(const char*const argv0, const int ret)
{
        std::cerr << "Usage: " << argv0 << 1+R"(
 [-r referenceSphereRadiusInKM] [-m metersPerUnitValue] -i inputFile.IMG -o outputFileName
Options:
 --help                 Show this help message and quit
 -r,--radius NUM        Radius of the reference sphere in km for the output
 -m,--m-per-unit NUM    Unit of altitude to use in the output
 -i FILE                Input file from SLDEM2015, must be accompanied by the corresponding .LBL file
 -o FILE                Output image path
)";
    return ret;
}


int main(int argc, char** argv)
try
{
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    double outputReferenceRadius = -1;
    double outputMetersPerUnit = -1;
    QString inFileName;
    QString outFileName;

    if(argc == 1)
        return usage(argv[0],1);

    for(int n = 1; n < argc; ++n)
    {
#define REQUIRE_PARAM() do{                                                 \
            if(argc+1 == n)                                                 \
            {                                                               \
                std::cerr << "Option " << arg << " requires parameter\n";   \
                return usage(argv[0], 1);                                   \
            }}while(false)

        const std::string arg = argv[n];
        if(arg == "--help")
        {
            return usage(argv[0], 0);
        }
        else if(arg == "-m" || arg == "--m-per-unit")
        {
            REQUIRE_PARAM();
            outputMetersPerUnit = std::stod(argv[++n]);
        }
        else if(arg == "-r" || arg == "--radius")
        {
            REQUIRE_PARAM();
            outputReferenceRadius = std::stod(argv[++n]);
        }
        else if(arg == "-i")
        {
            REQUIRE_PARAM();
            inFileName = argv[++n];
        }
        else if(arg == "-o")
        {
            REQUIRE_PARAM();
            outFileName = argv[++n];
        }
        else
        {
            std::cerr << "Unknown option " << arg << "\n";
            return usage(argv[0], 1);
        }
    }
#undef REQUIRE_PARAM

    if(inFileName.isEmpty() || outFileName.isEmpty())
    {
        std::cerr << (inFileName.isEmpty() ? "Input" : "Output") << " file was not specified\n";
        return usage(argv[0],1);
    }

    if(!inFileName.toUpper().endsWith(".IMG"))
    {
        std::cerr << "Input file name must end with .IMG\n";
        return 1;
    }
    const auto labelFileName = inFileName.endsWith(".IMG") ? QString(inFileName).replace(QRegularExpression("\\.IMG$"), ".LBL")
                                                           : QString(inFileName).replace(QRegularExpression("\\.img$"), ".lbl");

    QFile file(labelFileName);
    if(!file.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open " << inFileName.toStdString() << " for reading: "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }
    const auto label = file.readAll();
    if(label.isEmpty())
    {
        std::cerr << "Failed to read label file " << labelFileName.toStdString() << ": "
                  << file.errorString().toStdString() << "\n";
        return 1;
    }

    const auto width  = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();
    const auto height = QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(label).captured(1).toUInt();

    const auto radiusA = QRegularExpression(R"(\n\s*A_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();
    const auto radiusB = QRegularExpression(R"(\n\s*B_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();
    const auto radiusC = QRegularExpression(R"(\n\s*C_AXIS_RADIUS\s*=\s*([0-9.]+)\s*<(?:KM|km)>\s*\n)").match(label).captured(1).toDouble();

    const auto sampleType = QRegularExpression(R"(\n\s*\bSAMPLE_TYPE\s*=\s*([^\s]+)\s*\n)").match(label).captured(1);
    const auto sampleBits = QRegularExpression(R"(\n\s*\bSAMPLE_BITS\s*=\s*([^\s]+)\s*\n)").match(label).captured(1);
    const auto sampleOffset = QRegularExpression(R"(\n\s*\bOFFSET\s*=\s*([0-9.]+)\s*\n)").match(label).captured(1).toDouble();
    const auto sampleScalingFactor = QRegularExpression(R"(\bSCALING_FACTOR\s*=\s*([0-9.]+)\s*\n)").match(label).captured(1).toDouble();
    const auto sampleUnit = QRegularExpression(R"(\n\s*\bUNIT\s*=\s*([^\s]+)\b\s*\n)").match(label).captured(1);

    enum class Type
    {
        Float32,
        Int16,
    } inputType;

    if(sampleType == "PC_REAL" && sampleBits == "32")
    {
        inputType = Type::Float32;
    }
    else if(sampleType == "LSB_INTEGER" && sampleBits == "16")
    {
        inputType = Type::Int16;
    }
    else
    {
        throw std::runtime_error(QString("Sample type %1 with %2 bits per sample isn't supported. Only (little-endian) float32 and int16 are.")
                                    .arg(sampleType).arg(sampleBits).toStdString());
    }

    if(radiusA <= 0 || radiusA != radiusB || radiusB != radiusC)
    {
        throw std::runtime_error(QString("Bad reference sphere radii found. Expected three equal positive radii, but found %1, %2, %3")
                                    .arg(radiusA).arg(radiusB).arg(radiusC).toStdString());
    }

    if(width==0 || height==0)
    {
        std::cerr << "Failed to find input image dimensions\n";
        return 1;
    }

    if(sampleUnit != "KILOMETER" && sampleUnit != "METER")
        throw std::runtime_error(("Unexpected sample unit: "+sampleUnit).toStdString());

    const bool samplesAreInKM = sampleUnit == "KILOMETER";
    const auto sampleOffsetInKM = sampleOffset * (samplesAreInKM ? 1 : 1e-3);

    const auto inputMetersPerUnit = sampleScalingFactor * (samplesAreInKM ? 1000 : 1);
    if(outputMetersPerUnit <= 0)
        outputMetersPerUnit = 1;

    const auto inputReferenceRadius = sampleOffset;
    if(outputReferenceRadius <= 0)
        outputReferenceRadius = inputReferenceRadius * (samplesAreInKM ? 1 : 1000);

    if(std::abs(sampleOffsetInKM - inputReferenceRadius) > 1e-4)
    {
        throw std::runtime_error(QString("Unexpected sample offset: expected %1, got %2.")
                                    .arg(inputReferenceRadius).arg(sampleOffsetInKM).toStdString());
    }

    QFile in(inFileName);
    if(!in.open(QFile::ReadOnly))
    {
        std::cerr << "Failed to open input file for reading: "
                  << in.errorString().toStdString() << "\n";
        return 1;
    }

    QImage output;
    output = QImage(width, height, QImage::Format_Grayscale16);
    output.fill(QColor(0,0,0));
    const auto outData = reinterpret_cast<uint16_t*>(output.bits());
    const ssize_t numPixels = width * height;

    switch(inputType)
    {
    case Type::Float32:
    {
        const float* data = nullptr;
        const qint64 sizeToRead = width * height * sizeof data[0];
        data = reinterpret_cast<decltype(data)>(in.map(0, sizeToRead));
        for(ssize_t n = 0; n < numPixels; ++n)
        {
            const double value = data[n];
            const auto recenteredValue = value + inputReferenceRadius - outputReferenceRadius;
            const auto outputValue = recenteredValue * inputMetersPerUnit / outputMetersPerUnit;
            outData[n] = outputValue;
        }
        break;
    }
    case Type::Int16:
    {
        const int16_t* data = nullptr;
        const qint64 sizeToRead = width * height * sizeof data[0];
        data = reinterpret_cast<decltype(data)>(in.map(0, sizeToRead));
        for(ssize_t n = 0; n < numPixels; ++n)
        {
            const double value = data[n];
            const auto recenteredValue = value + inputReferenceRadius - outputReferenceRadius;
            const auto outputValue = recenteredValue * inputMetersPerUnit / outputMetersPerUnit;
            outData[n] = outputValue;
        }
        break;
    }
    }

    std::cerr << "Saving output file... ";
    if(!output.save(outFileName))
    {
        std::cerr << "Failed to save output file\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}

