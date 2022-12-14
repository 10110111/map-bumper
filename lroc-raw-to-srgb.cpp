#include <iostream>

#include <QFile>
#include <QImage>
#include <QByteArray>
#include <QRegularExpression>

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    return step(0.0031308,c)*(1.055*std::pow(c, 1/2.4)-0.055)+step(-0.0031308,-c)*12.92*c;
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

    // These values were computed from CIE 1931 functions using triangular
    // indicator functions (result should be equivalent to using linearly
    // interpolated spectra).
    constexpr double rgbsPerWL[][3] = {{ 0.0859932073132241,-0.0742345015484867 , 0.522206963569252 },
                                       {-9.8696160185158   , 14.3897262411038   , 81.0646570472351  },
                                       { 3.56460501999052  , 84.6557094421366   , 18.5537216323101  },
                                       { 83.5638258393946  , 6.38930071626876   ,-2.47243752916485  },
                                       { 45.0832972872074  ,-3.37526405235609   ,-0.491743508354297 },
                                       { 5.39818758379944  ,-0.559432795037625  ,-0.0396279820939476}};

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
    std::cerr << "Saving output file with dimensions " << width << " ?? " << height << "... ";
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
