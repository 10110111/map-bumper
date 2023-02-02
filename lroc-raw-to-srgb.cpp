#include <iostream>

#include <QFile>
#include <QImage>
#include <QByteArray>
#include <QRegularExpression>

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    const auto s = step(0.0031308,c);
    return s  * (1.055*std::pow(c, 1/2.4)-0.055) +
        (1-s) *  12.92*c;
}

int main(int argc, char** argv)
try
{
    if(argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " inputDir sector outputFile [valueScale]\n";
        std::cerr << "\"Sector\" is the geographic position part of the filename e.g. \"E350N0450\".\n";
        return 1;
    }

    const std::string inDir  = argv[1];
    const std::string sector = argv[2];
    const std::string outFileName = argv[3];
    const double valueScale = (argc==5 ? std::stod(argv[4]) : 1.) / 105.;

    constexpr int wavelengths[] = {360,415,566,604,643,689};
    constexpr size_t numWLs = std::size(wavelengths);
    std::vector<std::vector<float>> dataPerWL(numWLs);
    ssize_t width = -1, height = -1, rowStride = -1;
    for(size_t wlN=0; wlN<std::size(wavelengths); ++wlN)
    {
        const auto filename = inDir + "/WAC_HAPKE_" + std::to_string(wavelengths[wlN])
                                    + "NM_" + sector + ".IMG";
        QFile file(filename.c_str());
        if(!file.open(QFile::ReadOnly))
        {
            std::cerr << "Failed to open " << filename << " for reading: "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
        QByteArray header(27360, '\0');
        if(file.read(header.data(), header.size()) != header.size())
        {
            std::cerr << "Failed to read header of " << filename << ": "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
        const auto currWidth = QRegularExpression("SAMPLE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
        const auto currHeight= QRegularExpression("LINE_LAST_PIXEL\\s*=\\s*([0-9]+)\\b").match(header).captured(1).toUInt();
        if(currWidth==0 || currHeight==0)
        {
            std::cerr << "Failed to find image dimensions\n";
            return 1;
        }

        const auto currRowStride = currWidth;
        if(width<0)
            width = currWidth;
        if(height<0)
            height = currHeight;
        if(rowStride<0)
            rowStride = currRowStride;
        if(currWidth != width || currHeight != height || currRowStride != rowStride)
        {
            std::cerr << "Image dimensions or layout mismatch at " << filename << "\n";
            return 1;
        }

        dataPerWL[wlN].resize(rowStride*height);
        const qint64 sizeToRead = dataPerWL[wlN].size() * sizeof dataPerWL[wlN][0];
        if(file.read(reinterpret_cast<char*>(dataPerWL[wlN].data()), sizeToRead) != sizeToRead)
        {
            std::cerr << "Failed to read " << filename << ": "
                      << file.errorString().toStdString() << "\n";
            return 1;
        }
    }

    // These values were computed from CIE 1931 functions using triangular indicator functions
    // multiplied with solar irradiance at near-Earth orbit (result should be equivalent to using
    // linearly interpolated spectra of lunar albedo).
    constexpr double rgbsPerWL[][3] = {{0.123544388739659, -0.106861878682004, 0.752259001205208 },
                                       {-19.8366079265001,  27.8275045114262,  151.197577413019  },
                                       { 4.43497668468578,  157.373608612933,  37.079183527029   },
                                       { 145.775550325006,  11.7192895881426, -4.38388154872636  },
                                       { 73.8247972448887, -5.47353897213371, -0.811862808601448 },
                                       { 8.24022224911719, -0.85351964335178, -0.0605466729302762}};

    static_assert(std::size(rgbsPerWL) == numWLs);
    std::vector<double> out(rowStride*height*3);
    for(size_t n = 0; n < dataPerWL[0].size(); ++n)
    {
        double r = 0, g = 0, b = 0;
        bool good = true;
        for(size_t wlN = 0; wlN < numWLs; ++wlN)
        {
            const auto value = dataPerWL[wlN][n];
            if(value < -1e38 || !std::isfinite(value))
            {
                good = false;
                break;
            }
            r += value * rgbsPerWL[wlN][0];
            g += value * rgbsPerWL[wlN][1];
            b += value * rgbsPerWL[wlN][2];
        }
        if(good)
        {
            out[3*n+0] = r;
            out[3*n+1] = g;
            out[3*n+2] = b;
        }
        else
        {
            // Mark with magenta
            out[3*n+0] = 1/valueScale;
            out[3*n+1] = 0;
            out[3*n+2] = 1/valueScale;
        }
    }

    std::vector<uint8_t> outImgData(rowStride*height*3);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(out[n] * valueScale, 0., 1.));
    const QImage outImg(outImgData.data(), width, height, rowStride*3, QImage::Format_RGB888);
    std::cerr << "Saving output file with dimensions " << width << " × " << height << "... ";
    if(!outImg.save(outFileName.c_str()))
    {
        std::cerr << "Failed to save output file " << outFileName << "\n";
        return 1;
    }
    std::cerr << "done\n";
}
catch(std::exception const& ex)
{
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
}
