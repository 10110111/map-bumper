#include <iostream>

#include <QImage>

template<typename T> T step(T edge, T x) { return x<edge ? 0 : 1; }
double sRGBTransferFunction(const double c)
{
    return step(0.0031308,c)*(1.055*std::pow(c, 1/2.4)-0.055)+step(-0.0031308,-c)*12.92*c;
}

int main(int argc, char** argv)
try
{
    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " inputDir sector outputFile\n";
        std::cerr << "\"Sector\" is the geographic position part of the filename e.g. \"E350N0450\".\n";
        return 1;
    }

    const std::string inDir  = argv[1];
    const std::string sector = argv[2];
    const std::string outFileName = argv[3];

    constexpr int wavelengths[] = {360,415,566,604,643,689};
    std::vector<std::vector<uint8_t>> dataPerWL(std::size(wavelengths));
    ssize_t width = -1, height = -1, rowStride = -1;
    for(size_t n=0; n<std::size(wavelengths); ++n)
    {
        const auto filename = inDir + "/WAC_HAPKE_" + std::to_string(wavelengths[n]) + "NM_" + sector + ".tiff";
        QImage in(filename.c_str());
        if(in.isNull())
        {
            std::cerr << "Failed to read image from " << filename << "\n";
            return 1;
        }
        if(!in.isGrayscale())
        {
            std::cerr << "Input image " << filename << " is not grayscale\n";
            return 1;
        }
        in = in.convertToFormat(QImage::Format_Grayscale8);

        const auto currWidth  = in.width();
        const auto currHeight = in.height();
        const auto currRowStride = in.bytesPerLine();
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

        dataPerWL[n].resize(rowStride*height);
        std::copy_n(in.bits(), dataPerWL[n].size(), dataPerWL[n].data());
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

    static_assert(std::size(rgbsPerWL) == std::size(wavelengths));
    std::vector<double> out(rowStride*height*3);
    for(size_t n = 0; n < dataPerWL[0].size(); ++n)
    {
        double r = 0, g = 0, b = 0;
        for(size_t i = 0; i < dataPerWL.size(); ++i)
        {
            r += dataPerWL[i][n] * rgbsPerWL[i][0];
            g += dataPerWL[i][n] * rgbsPerWL[i][1];
            b += dataPerWL[i][n] * rgbsPerWL[i][2];
        }
        out[3*n+0] = r;
        out[3*n+1] = g;
        out[3*n+2] = b;
    }
    const auto max = *std::max_element(out.begin(), out.end());

    std::vector<uint8_t> outImgData(rowStride*height*3);
    for(size_t n = 0; n < out.size(); ++n)
        outImgData[n] = 255.*sRGBTransferFunction(std::clamp(out[n] / max, 0., 1.));
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
