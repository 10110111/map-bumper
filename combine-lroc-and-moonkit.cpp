#include <iostream>
#include <QImage>

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
    // Override the default limit of 128MiB, but let the user override our choice too
    setenv("QT_IMAGEIO_MAXALLOC","8192",false);

    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " moonkit.png lroc.png output.png\n";
        return 1;
    }

    const QString moonkitFileName  = argv[1];
    const QString lrocFileName  = argv[2];
    const QString outFileName = argv[3];

    std::cerr << "Loading LROC image...\n";
    const auto lroc = QImage(lrocFileName).convertToFormat(QImage::Format_RGB888);
    if(lroc.isNull())
    {
        std::cerr << "Failed to read LROC image\n";
        return 1;
    }

    std::cerr << "Loading Moon Kit image...\n";
    auto moonkit = QImage(moonkitFileName).convertToFormat(QImage::Format_RGB888);
    if(moonkit.isNull())
    {
        std::cerr << "Failed to read moonkit image\n";
        return 1;
    }

    std::cerr << "Combining...\n";
    const auto lrocShift = lroc.height() % 2; // skip the first, useless, line if present
    const int W = moonkit.width(), H = moonkit.height();
    const int W3 = W*3;
    const int moonkitRowStride = moonkit.bytesPerLine();
    const int lrocRowStride = lroc.bytesPerLine();
    uchar* data = moonkit.bits();
    const uchar* lrocData = lroc.bits() + lrocShift * lrocRowStride;
    const auto polesLowerBorder = (90-70)/180.*moonkit.height();
    const auto polesUpperBorder = polesLowerBorder + lroc.height() - lrocShift;

    for(int j = 0; j < H; ++j)
    {
        if(j < polesLowerBorder)
        {
            // North pole
            for(int i = 0; i < W3; i += 3)
            {
                constexpr double shiftR = 0.75182503818554 , coefR = -0.0523835897824206;
                constexpr double shiftG = 0.596639979156568, coefG = -0.0443378282579531;
                constexpr double shiftB = 0.4507008510799  , coefB = -0.0338486025950246;
                constexpr double power = -4.4;
                const auto u = data[i+1] / 255.;
                const auto v = std::pow(u, power);
                const auto r = std::clamp(shiftR + coefR * v, 0., 1.);
                const auto g = std::clamp(shiftG + coefG * v, 0., 1.);
                const auto b = std::clamp(shiftB + coefB * v, 0., 1.);
                data[i + 0] = sRGBTransferFunction(r) * 255;
                data[i + 1] = sRGBTransferFunction(g) * 255;
                data[i + 2] = sRGBTransferFunction(b) * 255;
            }
        }
        else if(j >= polesUpperBorder)
        {
            // South pole
            for(int i = 0, k = 0; i < W3; i += 3, ++k)
            {
                constexpr double shift = 0.603591562174297, coef = -0.0416058207644435, power = -4.4;
                const double valMK = data[i+1] / 255.;
                const double Y = std::max(0.24, shift + coef * std::pow(valMK, power));
                const double r0 = 1.26443109577546, g0 = 0.950351679150597, b0 = 0.712763115940595;
                const double r = r0 * Y, g = g0 * Y, b = b0 * Y;
                data[i + 0] = sRGBTransferFunction(std::clamp(r, 0., 1.)) * 255;
                data[i + 1] = sRGBTransferFunction(std::clamp(g, 0., 1.)) * 255;
                data[i + 2] = sRGBTransferFunction(std::clamp(b, 0., 1.)) * 255;
            }
        }
        else
        {
            std::memcpy(data, lrocData, W3);
            lrocData += lrocRowStride;
        }
        data += moonkitRowStride;
    }

    std::cerr << "Saving output file... ";
    if(!moonkit.save(outFileName))
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
