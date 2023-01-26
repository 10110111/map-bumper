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
    for(int n = 0; n < numBytes; ++n)
    {
        const auto u = data[n] / 255.;
        const auto a = std::clamp(0.6*u+0.4, 0., 1.);
        const auto b = std::pow(a, 2.8);
        data[n] = sRGBTransferFunction(b) * 255;
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
