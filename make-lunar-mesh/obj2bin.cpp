#include <vector>
#include <fstream>
#include <stdint.h>
#include <iostream>

constexpr double AU = 149597870.691; // km

struct Model
{
    std::vector<float> vertices;
    std::vector<uint32_t> indices;
};

Model readOBJ(const char*const filename)
{
    std::ifstream file(filename);
    if(!file)
        throw std::runtime_error("Failed to open input file");

    Model model;
    std::string line;
    while(std::getline(file, line))
    {
        if(line.empty()) continue;
        switch(line[0])
        {
        case 'v':
        {
            if(line.size() >= 2 && line[1] != ' ')
                continue;
            float x, y, z;
            if(sscanf(line.c_str(), "v %g %g %g", &x, &y, &z) != 3)
                throw std::runtime_error("Bad vertex line: \""+line+"\"");
            model.vertices.push_back(x / AU);
            model.vertices.push_back(y / AU);
            model.vertices.push_back(z / AU);
            break;
        }
        case 'f':
        {
            int i1, i2, i3;
            if(sscanf(line.c_str(), "f %d %d %d", &i1, &i2, &i3) != 3)
                throw std::runtime_error("Bad or unsupported face line: \""+line+"\"");
            model.indices.push_back(i1);
            model.indices.push_back(i2);
            model.indices.push_back(i3);
            break;
        }
        }
    }

    return model;
}

int main(int argc, char** argv)
try
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " filename.obj filename.bin\n";
        return 1;
    }
    const char*const objFileName = argv[1];
    const char*const binFileName = argv[2];
    const auto model = readOBJ(objFileName);

    std::cerr << "OBJ vertices: " << model.vertices.size() / 3 << "\n";
    std::cerr << "OBJ faces: " << model.indices.size() / 3 << "\n";

    std::ofstream out(binFileName, std::ios_base::binary);
    const uint32_t vertexCount = model.vertices.size() / 3;
    const uint32_t indexCount = model.indices.size();
    out.write(reinterpret_cast<const char*>(&vertexCount), sizeof vertexCount);
    out.write(reinterpret_cast<const char*>(&indexCount), sizeof indexCount);
    out.write(reinterpret_cast<const char*>(model.vertices.data()),
              model.vertices.size()*sizeof model.vertices[0]);
    out.write(reinterpret_cast<const char*>(model.indices.data()),
              model.indices.size()*sizeof model.indices[0]);
    out.close();
    if(!out)
    {
        std::cerr << "Failed to write output file\n";
        return 1;
    }
}
catch(std::exception const& ex)
{
    std::cerr << ex.what() << "\n";
    return 1;
}
