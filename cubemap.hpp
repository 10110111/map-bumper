#pragma once
#include <stdint.h>

static constexpr int MAX_ALT_TILES_PER_CUBE_SIDE = 4;
static constexpr int CUBEMAP_FORMAT_VERSION = 3;
struct CubeMapFileHeader
{
    uint32_t formatVersion;
    uint32_t cubeMapSide;
    int16_t maxAltitudes[MAX_ALT_TILES_PER_CUBE_SIDE][MAX_ALT_TILES_PER_CUBE_SIDE];
};
