#include "maths.h"
#include <stdlib.h>
#include <stdint.h>

float RandomFloat01(unsigned int *seed0, unsigned int *seed1)
{
    *seed0 = 36969 * ((*seed0) & 65535) + ((*seed0) >> 16);
    *seed1 = 18000 * ((*seed1) & 65535) + ((*seed1) >> 16);

    unsigned int ires = ((*seed0) << 16) + (*seed1);

    /* Convert to float */
    union {
        float f;
        unsigned int ui;
    } res;
    res.ui = (ires & 0x007fffff) | 0x40000000;

    return (res.f - 2.f) / 2.f;
}

float3 RandomInUnitDisk(unsigned int *seed0, unsigned int *seed1)
{
    float px, py;
    do {
        px = 2.0f*RandomFloat01(seed0, seed1) - 1.0f;
        py = 2.0f*RandomFloat01(seed0, seed1) - 1.0f;
    } while (px*px + py*py >= 1.0f);

    return float3(px, py, 0.0f);
}

float3 RandomUnitVector(unsigned int *seed0, unsigned int *seed1)
{
    float px, py, pz;

    do {
        px = 2.0f*RandomFloat01(seed0, seed1) - 1.0f;
        py = 2.0f*RandomFloat01(seed0, seed1) - 1.0f;
        pz = 2.0f*RandomFloat01(seed0, seed1) - 1.0f;
    } while (px*px + py*py + pz*pz >= 1.0f);

    return float3(px, py, pz);
}
