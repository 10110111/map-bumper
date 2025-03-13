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
    setenv("QT_IMAGEIO_MAXALLOC","4096",false);

    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " input output\n";
        return 1;
    }

    const QString inFileName  = argv[1];
    const QString outFileName = argv[2];

    auto img = QImage(inFileName).convertToFormat(QImage::Format_RGB888);
    if(img.isNull())
    {
        std::cerr << "Failed to read input image\n";
        return 1;
    }
    uchar*const data = img.bits();
    const auto numBytes = img.height() * img.bytesPerLine();
    for(int n = 0; n < numBytes; n += 3)
    {
        constexpr double shiftR = 0.75182503818554 , coefR = -0.0523835897824206;
        constexpr double shiftG = 0.596639979156568, coefG = -0.0443378282579531;
        constexpr double shiftB = 0.4507008510799  , coefB = -0.0338486025950246;
        constexpr double power = -4.4;
        const auto u = data[n] / 255.;
        const auto v = std::pow(u, power);
        const auto r = std::clamp(shiftR + coefR * v, 0., 1.);
        const auto g = std::clamp(shiftG + coefG * v, 0., 1.);
        const auto b = std::clamp(shiftB + coefB * v, 0., 1.);
        data[n + 0] = sRGBTransferFunction(r) * 255;
        data[n + 1] = sRGBTransferFunction(g) * 255;
        data[n + 2] = sRGBTransferFunction(b) * 255;
    }
    std::cerr << "Saving output file... ";
    if(!img.save(outFileName))
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
