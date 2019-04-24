// main program entry point and the actual raytracing bits

#include "maths.h"
#include "scene.h"


// Include external libraries:
// - PNG writing
#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"
// - time measurement
#define SOKOL_IMPL
#include "external/sokol_time.h"
// - OBJ file loading
#include "external/objparser.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _OPENCL
#include "external/clbase.h"
#include "external/CLEW/clew.h"
#endif

#include <vector>
#include <numeric>

// --------------------------------------------------------------------------
// "ray/path tracing" bits

// general minimum/maximum distances for rays (from "very close to surface but not exacttly on it"
// to "ten million units")
const float kMinT = 0.001f;
const float kMaxT = 1.0e7f;
// maximum raytracing recursion depth, i.e. number of light bounces
const int kMaxDepth = 10;

// we have one hardcoded directional light, with this direction and color
static const float3 kLightDir = normalize(float3(-0.7f,1.0f,0.5f));
static const float3 kLightColor = float3(0.7f,0.6f,0.5f);


// when a ray "r" has just hit a surface at point "hit", decide what to do about it:
// in our very simple case, we assume the surface is perfectly diffuse, so we'll return:
// - surface albedo ("color") in "attenuation"
// - new random ray for the next light bounce in "scattered"
// - illumination from the directional light in "outLightE"
static bool Scatter(const Ray& r, const Hit& hit, float3& attenuation, Ray& scattered, float3& outLightE, unsigned int *seed0, unsigned int *seed1, int& inoutRayCount)
{
    outLightE = float3(0,0,0);

    // model a perfectly diffuse material:
    
    // random point on unit sphere that is tangent to the hit point
    float3 target = hit.pos + hit.normal + RandomUnitVector(seed0, seed1);
    scattered = Ray(hit.pos, normalize(target - hit.pos));
    
    // make color slightly based on surface normals
    float3 albedo = hit.normal * 0.0f + float3(0.7f,0.7f,0.7f);
    attenuation = albedo;
    
    // explicit directional light by shooting a shadow ray
    ++inoutRayCount;

    Hit lightHit;
    int id = HitScene(Ray(hit.pos, kLightDir), kMinT, kMaxT, lightHit);
    if (id == -1)
    {
        // ray towards the light did not hit anything in the scene, so
        // that means we are not in shadow: compute illumination from it
        float3 rdir = r.dir;
        AssertUnit(rdir);
        float3 nl = dot(hit.normal, rdir) < 0 ? hit.normal : -hit.normal;
        outLightE += albedo * kLightColor * (fmax(0.0f, dot(kLightDir, nl)));
    }

    return true;
}


// trace a ray into the scene, and return the normal color for it (debug)
static float3 TraceNormal(const Ray& r, int depth, int& inoutRayCount)
{
    (void)depth;

    ++inoutRayCount;

    Hit hit;
    int id = HitScene(r, kMinT, kMaxT, hit);
    if (id != -1)
    {
        // ray hits something in the scene: return world normal
        return float3(std::abs(hit.normal.getX()), std::abs(hit.normal.getY()), std::abs(hit.normal.getZ()));
    }
    else
    {
        // ray does not hit anything: return black
        return float3(0, 0, 0);
    }
}


// trace a ray into the scene, and return the final color for it
static float3 Trace(const Ray& r, int depth, unsigned int *seed0, unsigned int *seed1, int& inoutRayCount)
{
    ++inoutRayCount;

    Hit hit;
    int id = HitScene(r, kMinT, kMaxT, hit);
    if (id != -1)
    {
        // ray hits something in the scene
        Ray scattered;
        float3 attenuation;
        float3 lightE;
        if (depth < kMaxDepth && Scatter(r, hit, attenuation, scattered, lightE, seed0, seed1, inoutRayCount))
        {
            // we got a new ray bounced from the surface; recursively trace it
            return lightE + attenuation * Trace(scattered, depth+1, seed0, seed1, inoutRayCount);
        }
        else
        {
            // reached recursion limit, or surface fully absorbed the ray: return black
            return float3(0,0,0);
        }
    }
    else
    {
        // ray does not hit anything: return illumination from the sky (just a simple gradient really)
        float3 unitDir = r.dir;
        float t = 0.5f*(unitDir.getY() + 1.0f);
        return ((1.0f - t)*float3(1.0f, 1.0f, 1.0f) + t * float3(0.5f, 0.7f, 1.0f)) * 0.5f;
    }
}


// trace a ray into the scene, and return the final color for it
static float3 TraceIterative(const Ray& r, int depth, unsigned int *seed0, unsigned int *seed1, int& inoutRayCount)
{
    (void)depth;

    Ray currRay = r;
    float3 currAttenuation(1, 1, 1);
    float3 currLightE(0, 0, 0);

    for (int currDepth = 0; currDepth < kMaxDepth; ++currDepth)
    {
        ++inoutRayCount;

        Hit hit;
        int id = HitScene(currRay, kMinT, kMaxT, hit);
        if (id != -1)
        {
            // ray hits something in the scene
            Ray scattered;
            float3 attenuation;
            float3 lightE;
            if (Scatter(currRay, hit, attenuation, scattered, lightE, seed0, seed1, inoutRayCount))
            {
                currLightE += lightE * currAttenuation;
                currAttenuation *= attenuation;

                currRay = scattered;
            }
            else
            {
                // reached recursion limit, or surface fully absorbed the ray: return black
                return float3(0, 0, 0);
            }
        }
        else
        {
            // ray does not hit anything: return illumination from the sky (just a simple gradient really)
            float3 unitDir = currRay.dir;
            float t = 0.5f*(unitDir.getY() + 1.0f);
            float3 c = ((1.0f - t)*float3(1.0f, 1.0f, 1.0f) + t * float3(0.5f, 0.7f, 1.0f)) * 0.5f;
            return currLightE + currAttenuation * c;
        }
    }

    return float3(0, 0, 0); // exceeded recursion
}


// load scene from an .OBJ file
static bool LoadScene(const char* dataFile, float3& outBoundsMin, float3& outBoundsMax, std::vector<float>& trisFloat4)
{
    ObjFile objFile;
    if (!objParseFile(objFile, dataFile))
    {
        printf("ERROR: failed to load .obj file\n");
        return false;
    }
    outBoundsMin = float3(+1.0e6f, +1.0e6f, +1.0e6f);
    outBoundsMax = float3(-1.0e6f, -1.0e6f, -1.0e6f);

    int objTriCount = int(objFile.f_size / 9);
    Triangle* tris = new Triangle[objTriCount + 2]; // will add two triangles for the "floor"
    for (int i = 0; i < objTriCount; ++i)
    {
        int idx0 = objFile.f[i * 9 + 0] * 3;
        int idx1 = objFile.f[i * 9 + 3] * 3;
        int idx2 = objFile.f[i * 9 + 6] * 3;
        float3 v0 = float3(objFile.v[idx0 + 0], objFile.v[idx0 + 1], objFile.v[idx0 + 2]);
        float3 v1 = float3(objFile.v[idx1 + 0], objFile.v[idx1 + 1], objFile.v[idx1 + 2]);
        float3 v2 = float3(objFile.v[idx2 + 0], objFile.v[idx2 + 1], objFile.v[idx2 + 2]);
        tris[i].v0 = v0;
        tris[i].v1 = v1;
        tris[i].v2 = v2;
        outBoundsMin = min(outBoundsMin, v0); outBoundsMax = max(outBoundsMax, v0);
        outBoundsMin = min(outBoundsMin, v1); outBoundsMax = max(outBoundsMax, v1);
        outBoundsMin = min(outBoundsMin, v2); outBoundsMax = max(outBoundsMax, v2);

        {
            // ocl memory helpers

            trisFloat4.emplace_back(tris[i].v0.getX());
            trisFloat4.emplace_back(tris[i].v0.getY());
            trisFloat4.emplace_back(tris[i].v0.getZ());
            trisFloat4.emplace_back(0.0f);

            trisFloat4.emplace_back(tris[i].v1.getX());
            trisFloat4.emplace_back(tris[i].v1.getY());
            trisFloat4.emplace_back(tris[i].v1.getZ());
            trisFloat4.emplace_back(0.0f);

            trisFloat4.emplace_back(tris[i].v2.getX());
            trisFloat4.emplace_back(tris[i].v2.getY());
            trisFloat4.emplace_back(tris[i].v2.getZ());
            trisFloat4.emplace_back(0.0f);
        }
    }

    // add two triangles that are right "under the scene" and covering larger area than the scene
    // itself, to serve as a "floor"
    float3 size = outBoundsMax - outBoundsMin;
    float3 extra = size * 0.7f;
    tris[objTriCount+0].v0 = float3(outBoundsMin.x-extra.x, outBoundsMin.y, outBoundsMin.z-extra.z);
    tris[objTriCount+0].v1 = float3(outBoundsMin.x-extra.x, outBoundsMin.y, outBoundsMax.z+extra.z);
    tris[objTriCount+0].v2 = float3(outBoundsMax.x+extra.x, outBoundsMin.y, outBoundsMin.z-extra.z);
    tris[objTriCount+1].v0 = float3(outBoundsMin.x-extra.x, outBoundsMin.y, outBoundsMax.z+extra.z);
    tris[objTriCount+1].v1 = float3(outBoundsMax.x+extra.x, outBoundsMin.y, outBoundsMax.z+extra.z);
    tris[objTriCount+1].v2 = float3(outBoundsMax.x+extra.x, outBoundsMin.y, outBoundsMin.z-extra.z);

    {
        // ocl memory helpers

        trisFloat4.emplace_back(tris[objTriCount+0].v0.getX());
        trisFloat4.emplace_back(tris[objTriCount+0].v0.getY());
        trisFloat4.emplace_back(tris[objTriCount+0].v0.getZ());
        trisFloat4.emplace_back(0.0f);

        trisFloat4.emplace_back(tris[objTriCount+0].v1.getX());
        trisFloat4.emplace_back(tris[objTriCount+0].v1.getY());
        trisFloat4.emplace_back(tris[objTriCount+0].v1.getZ());
        trisFloat4.emplace_back(0.0f);

        trisFloat4.emplace_back(tris[objTriCount+0].v2.getX());
        trisFloat4.emplace_back(tris[objTriCount+0].v2.getY());
        trisFloat4.emplace_back(tris[objTriCount+0].v2.getZ());
        trisFloat4.emplace_back(0.0f);

        trisFloat4.emplace_back(tris[objTriCount+1].v0.getX());
        trisFloat4.emplace_back(tris[objTriCount+1].v0.getY());
        trisFloat4.emplace_back(tris[objTriCount+1].v0.getZ());
        trisFloat4.emplace_back(0.0f);

        trisFloat4.emplace_back(tris[objTriCount+1].v1.getX());
        trisFloat4.emplace_back(tris[objTriCount+1].v1.getY());
        trisFloat4.emplace_back(tris[objTriCount+1].v1.getZ());
        trisFloat4.emplace_back(0.0f);

        trisFloat4.emplace_back(tris[objTriCount+1].v2.getX());
        trisFloat4.emplace_back(tris[objTriCount+1].v2.getY());
        trisFloat4.emplace_back(tris[objTriCount+1].v2.getZ());
        trisFloat4.emplace_back(0.0f);
    }

    uint64_t t0 = stm_now();
    InitializeScene(objTriCount + 2, tris);
    printf("Initialized scene '%s' (%i tris) in %.3fs\n", dataFile, objTriCount+2, stm_sec(stm_since(t0)));

    delete[] tris;
    return true;
}

struct TraceData
{
    int screenWidth, screenHeight, samplesPerPixel;
    std::vector<unsigned int> seed0;
    std::vector<unsigned int> seed1;
    uint8_t* image;
    const Camera* camera;
    int rayCount;
};

static void TraceImage(TraceData& data)
{
    uint8_t* image = data.image;
    float invWidth = 1.0f / data.screenWidth;
    float invHeight = 1.0f / data.screenHeight;

    std::vector<int> rayCountVec(data.screenHeight, 0);

    // go over the image: each pixel row
    #pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < data.screenHeight; ++y)
    {
        int rayCount = 0;

        // go over the image: each pixel in the row
        for (int x = 0; x < data.screenWidth; ++x)
        {
            uint32_t seed0 = data.seed0[x + y * data.screenWidth];
            uint32_t seed1 = data.seed1[x + y * data.screenWidth];

            float3 col(0, 0, 0);
            // we'll trace N slightly jittered rays for each pixel, to get anti-aliasing, loop over them here
            for (int s = 0; s < data.samplesPerPixel; ++s)
            {
                // get a ray from camera, and trace it
                float u = float(x + RandomFloat01(&seed0, &seed1)) * invWidth;
                float v = float(y + RandomFloat01(&seed0, &seed1)) * invHeight;
                Ray r = data.camera->GetRay(u, v, &seed0, &seed1);

                //col += Trace(r, 0, &seed0, &seed1, rayCount);
                //col += TraceIterative(r, 0, &seed0, &seed1, rayCount);
                col += TraceNormal(r, 0, rayCount);
            }
            col *= 1.0f / float(data.samplesPerPixel);

            // simplistic "gamma correction" by just taking a square root of the final color
            col.x = sqrtf(col.x);
            col.y = sqrtf(col.y);
            col.z = sqrtf(col.z);

            // our image is bytes in 0-255 range, turn our floats into them here and write into the image
            int index = x * 4 + data.screenWidth * 4 * y;
            image[index + 0] = uint8_t(saturate(col.x) * 255.0f);
            image[index + 1] = uint8_t(saturate(col.y) * 255.0f);
            image[index + 2] = uint8_t(saturate(col.z) * 255.0f);
            image[index + 3] = 255;
        }
        rayCountVec[y] = rayCount;
    }

    data.rayCount = std::accumulate(rayCountVec.begin(), rayCountVec.end(), 0);
}

#ifdef _OPENCL
// file path helper for OpenCL kernel file loading
bool findFullPath(const std::string& root, std::string& filePath)
{
    bool fileFound = false;
    const std::string resourcePath = root;

    filePath = resourcePath + filePath;
    for (unsigned int i = 0; i < 16; ++i)
    {
        std::ifstream file;
        file.open(filePath.c_str());
        if (file.is_open())
        {
            fileFound = true;
            break;
        }

        filePath = "../" + filePath;
    }

    return fileFound;
}
#endif

int main(int argc, const char** argv)
{
    std::string pngSuffix = "";

#ifdef _OPENMP
    pngSuffix = "_omp";

    int thread_num = omp_get_max_threads();
    printf("omp_get_max_threads: %d\n", thread_num);
#endif

#ifdef _OPENCL
    pngSuffix = "_ocl";

    int clewOK = initClew();
    if (clewOK != 0)
    {
        printf("ERROR: initClew() failed\n");
        return 1;
    }

    cl_context clContext = 0;
    cl_command_queue clQueue = 0;

    cl_program clProg = 0;
    cl_kernel clKernel = 0;

    cl_mem clMImage = 0;
    cl_mem clMTris = 0;
    cl_mem clMCamera = 0;
    cl_mem clMSeed0 = 0;
    cl_mem clMSeed1 = 0;

    OpenCLUtil cl;
    cl.init();
#endif

    // initialize timer
    stm_setup();

    // parse screen size command line arguments
    int screenWidth, screenHeight, samplesPerPixel;
    if (argc < 5)
    {
        printf("Usage: TrimeshTracer.exe [width] [height] [samplesPerPixel] [objFile]\n");
        return 1;
    }
    screenWidth = atoi(argv[1]);
    if (screenWidth < 1 || screenWidth > 10000)
    {
        printf("ERROR: invalid width argument '%s'\n", argv[1]);
        return 1;
    }
    screenHeight = atoi(argv[2]);
    if (screenHeight < 1 || screenHeight > 10000)
    {
        printf("ERROR: invalid height argument '%s'\n", argv[2]);
        return 1;
    }
    samplesPerPixel = atoi(argv[3]);
    if (samplesPerPixel < 1 || samplesPerPixel > 1024)
    {
        printf("ERROR: invalid samplesPerPixel argument '%s'\n", argv[3]);
        return 1;
    }

    // load model file and initialize the scene
    float3 sceneMin, sceneMax;
    std::vector<float> triOCLVec;
    if (!LoadScene(argv[4], sceneMin, sceneMax, triOCLVec))
        return 1;

    // place a camera: put it a bit outside scene bounds, looking at the center of it
    float3 sceneSize = sceneMax - sceneMin;
    float3 sceneCenter = (sceneMin + sceneMax) * 0.5f;
    float3 lookfrom = sceneCenter + sceneSize * float3(0.3f,0.6f,1.2f);
    if (strstr(argv[4], "sponza.obj") != nullptr) // sponza looks bad when viewed from outside; hardcode camera position
        lookfrom = float3(-5.96f, 4.08f, -1.22f);
    float3 lookat = sceneCenter + sceneSize * float3(0,-0.1f,0);
    float distToFocus = length(lookfrom - lookat);
    float aperture = 0.03f;
    auto camera = Camera(lookfrom, lookat, float3(0, 1, 0), 60, float(screenWidth) / float(screenHeight), aperture, distToFocus);

    std::vector<float> camOCLVec(8*4, 0.0f);
    {
        camOCLVec[ 0] = camera.origin.x;
        camOCLVec[ 1] = camera.origin.y;
        camOCLVec[ 2] = camera.origin.z;
        
        camOCLVec[ 4] = camera.lowerLeftCorner.x;
        camOCLVec[ 5] = camera.lowerLeftCorner.y;
        camOCLVec[ 6] = camera.lowerLeftCorner.z;

        camOCLVec[ 8] = camera.horizontal.x;
        camOCLVec[ 9] = camera.horizontal.y;
        camOCLVec[10] = camera.horizontal.z;

        camOCLVec[12] = camera.vertical.x;
        camOCLVec[13] = camera.vertical.y;
        camOCLVec[14] = camera.vertical.z;

        camOCLVec[16] = camera.u.x;
        camOCLVec[17] = camera.u.y;
        camOCLVec[18] = camera.u.z;

        camOCLVec[20] = camera.v.x;
        camOCLVec[21] = camera.v.y;
        camOCLVec[22] = camera.v.z;

        camOCLVec[24] = camera.w.x;
        camOCLVec[25] = camera.w.y;
        camOCLVec[26] = camera.w.z;

        camOCLVec[28] = camera.lensRadius;
        camOCLVec[29] = 0.0f;
        camOCLVec[30] = 0.0f;
        camOCLVec[31] = 0.0f;
    }

    std::vector<unsigned int> seed0OCLVec(screenWidth*screenHeight);
    std::vector<unsigned int> seed1OCLVec(screenWidth*screenHeight);
    {
        for (int i = 0; i < screenWidth*screenHeight; ++i)
        {
            seed0OCLVec[i] = rand() % RAND_MAX + 1;
            seed1OCLVec[i] = rand() % RAND_MAX + 1;
        }
    }

    // create RGBA image for the result
    uint8_t* image = new uint8_t[screenWidth * screenHeight * 4];

#ifdef _OPENCL
    {
        cl_uint selectedPlatformIndex;
        cl_uint selectedDeviceIndex;

        cl_platform_id selectedPlatformID;
        cl_device_id selectedDeviceID;

        // OpenCL context, device, queue init
        //cl.selectePlatformDevice_FirstCPU(selectedPlatformIndex, selectedDeviceIndex, selectedPlatformID, selectedDeviceID);
        cl.selectePlatformDevice_FirstGPU(selectedPlatformIndex, selectedDeviceIndex, selectedPlatformID, selectedDeviceID);

        if (selectedPlatformIndex == -1 || selectedDeviceIndex == -1)
        {
            printf("ERROR: could not get proper OpenCL device\n");
            return 1;
        }
        else
        {
            printf("OCL Platform: %s\n", cl.platforms[selectedPlatformIndex].platformName.c_str());
            printf("OCL Device: %s\n", cl.platforms[selectedPlatformIndex].devices[selectedDeviceIndex].deviceName.c_str());
        }
        
        std::vector<cl_context_properties> clContextProps = createContextProps(selectedPlatformID);

        cl_int clStatus = CL_SUCCESS;

        cl_int deviceCount = 1;
        clContext = clCreateContext(&clContextProps[0], deviceCount, &selectedDeviceID, NULL, NULL, &clStatus);
        CHECK_CL(clStatus);

        clQueue = clCreateCommandQueue(clContext, selectedDeviceID, 0, &clStatus);
        CHECK_CL(clStatus);

        // OpenCL kernel/program handling, memory management
        std::string rootStr = "data/";
        std::string filePath = "trimeshtracer.cl";

        bool fileFound = findFullPath(rootStr, filePath);
        if (fileFound == false)
        {
            printf("ERROR: could not find OpenCL kernel file\n");
            return 1;
        }

        std::ifstream myfile(filePath.c_str());
        std::string clProgStr(std::istreambuf_iterator<char>(myfile), (std::istreambuf_iterator<char>()));

        const char* clProgramStr = clProgStr.c_str();
        size_t clProgramSize = clProgStr.size();

        clProg = clCreateProgramWithSource(clContext, 1, &clProgramStr, &clProgramSize, &clStatus);
        CHECK_CL(clStatus);

        //clStatus = clBuildProgram(clProg_root, 0, NULL, NULL, NULL, NULL);
        clStatus = clBuildProgram(clProg, 0, NULL, "-cl-fast-relaxed-math", NULL, NULL);
        if (clStatus == CL_BUILD_PROGRAM_FAILURE)
        {
            size_t buildLogSize = 0;
            clStatus = clGetProgramBuildInfo(clProg, selectedDeviceID, CL_PROGRAM_BUILD_LOG, 0, NULL, &buildLogSize);
            CHECK_CL(clStatus);

            char* buildLog = new char[buildLogSize];
            clStatus = clGetProgramBuildInfo(clProg, selectedDeviceID, CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, NULL);
            CHECK_CL(clStatus);

            std::cout << buildLog << std::endl;
        }

        clKernel = clCreateKernel(clProg, "trace", &clStatus);
        CHECK_CL(clStatus);

        clMImage = clCreateBuffer(clContext, CL_MEM_READ_ONLY, 4 * screenWidth * screenHeight * sizeof(uint8_t), nullptr, &clStatus);
        CHECK_CL(clStatus);

        clMTris = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, triOCLVec.size() * sizeof(float), triOCLVec.data(), &clStatus);
        CHECK_CL(clStatus);

        clMCamera = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, camOCLVec.size() * sizeof(float), camOCLVec.data(), &clStatus);
        CHECK_CL(clStatus);

        clMSeed0 = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, screenWidth * screenHeight * sizeof(unsigned int), seed0OCLVec.data(), &clStatus);
        CHECK_CL(clStatus);

        clMSeed1 = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, screenWidth * screenHeight * sizeof(unsigned int), seed1OCLVec.data(), &clStatus);
        CHECK_CL(clStatus);
    }
#endif

    TraceData data;
    data.screenWidth = screenWidth;
    data.screenHeight = screenHeight;
    data.seed0 = seed0OCLVec;
    data.seed1 = seed1OCLVec;
    data.samplesPerPixel = samplesPerPixel;
    data.image = image;
    data.camera = &camera;
    data.rayCount = 0;

    // generate the image - run TraceImage
    uint64_t t0 = stm_now();

#ifndef _OPENCL
    TraceImage(data);
#else
    {
        cl_int clStatus = CL_SUCCESS;

        cl_int clSamplesPerPixel = samplesPerPixel;
        cl_int clWidth = screenWidth;
        cl_int clHeight = screenHeight;
        cl_float clInvWidth = 1.0f / screenWidth;
        cl_float clInvHeight = 1.0f / screenHeight;
        cl_int clTriangleNum = (cl_int)(triOCLVec.size()) / (3*4);

        size_t globalWS[3] = { (size_t)(screenWidth), (size_t)(screenHeight), 1 };
        size_t localWS[3] = { 8, 8, 1 };
        //size_t* localWS = NULL;

        clStatus |= clSetKernelArg(clKernel,  0, sizeof(cl_int), (void*)&clSamplesPerPixel);
        clStatus |= clSetKernelArg(clKernel,  1, sizeof(cl_int), (void*)&clWidth);
        clStatus |= clSetKernelArg(clKernel,  2, sizeof(cl_int), (void*)&clHeight);
        clStatus |= clSetKernelArg(clKernel,  3, sizeof(cl_float), (void*)&clInvWidth);
        clStatus |= clSetKernelArg(clKernel,  4, sizeof(cl_float), (void*)&clInvHeight);
        clStatus |= clSetKernelArg(clKernel,  5, sizeof(cl_int), (void*)&clTriangleNum);
        clStatus |= clSetKernelArg(clKernel,  6, sizeof(cl_mem), (void*)&clMTris);
        clStatus |= clSetKernelArg(clKernel,  7, sizeof(cl_mem), (void*)&clMSeed0);
        clStatus |= clSetKernelArg(clKernel,  8, sizeof(cl_mem), (void*)&clMSeed1);
        clStatus |= clSetKernelArg(clKernel,  9, sizeof(cl_mem), (void*)&clMCamera);
        clStatus |= clSetKernelArg(clKernel, 10, sizeof(cl_mem), (void*)&clMImage);
        CHECK_CL(clStatus);

        clStatus = clEnqueueNDRangeKernel(clQueue, clKernel, 2, NULL, globalWS, localWS, 0, NULL, NULL);
        CHECK_CL(clStatus);

        clStatus = clFinish(clQueue);
        CHECK_CL(clStatus);

        // skip first api init call
        t0 = stm_now();
        clStatus = clEnqueueNDRangeKernel(clQueue, clKernel, 2, NULL, globalWS, localWS, 0, NULL, NULL);
        CHECK_CL(clStatus);

        clStatus = clFinish(clQueue);
        CHECK_CL(clStatus);
    }
#endif    

    double dt = stm_sec(stm_since(t0));
    
#ifdef _OPENCL
    {
        cl_int clStatus = CL_SUCCESS;

        clStatus = clEnqueueReadBuffer(clQueue, clMImage, CL_TRUE, 0, screenWidth * screenHeight * 4 * sizeof(uint8_t), image, 0, nullptr, nullptr);
        CHECK_CL(clStatus);
    }
#endif

    printf("Rendered scene at %ix%i,%ispp in %.3f s\n", screenWidth, screenHeight, samplesPerPixel, dt);
    printf("- %.1f K Rays, %.1f K Rays/s\n", data.rayCount/1000.0, data.rayCount/1000.0/dt);

    // write resulting image as PNG
    stbi_flip_vertically_on_write(1);
    stbi_write_png((std::string("output") + pngSuffix + ".png").c_str(), screenWidth, screenHeight, 4, image, screenWidth*4);

    // cleanup and exit
    delete[] image;
    CleanupScene();

#ifdef _OPENCL
    {
        if (clMImage)
        {
            CHECK_CL(clReleaseMemObject(clMImage));
            clMImage = 0;
        }

        if (clMTris)
        {
            CHECK_CL(clReleaseMemObject(clMTris));
            clMTris = 0;
        }

        if (clMSeed0)
        {
            CHECK_CL(clReleaseMemObject(clMSeed0));
            clMSeed0 = 0;
        }

        if (clMSeed1)
        {
            CHECK_CL(clReleaseMemObject(clMSeed1));
            clMSeed1 = 0;
        }

        if (clMCamera)
        {
            CHECK_CL(clReleaseMemObject(clMCamera));
            clMCamera = 0;
        }

        if (clKernel)
        {
            CHECK_CL(clReleaseKernel(clKernel));
            clKernel = 0;
        }

        if (clProg)
        {
            CHECK_CL(clReleaseProgram(clProg));
            clProg = 0;
        }

        if (clQueue)
        {
            CHECK_CL(clReleaseCommandQueue(clQueue));
            clQueue = 0;
        }

        if (clContext)
        {
            CHECK_CL(clReleaseContext(clContext));
            clContext = 0;
        }
    }
#endif

    return 0;
}
