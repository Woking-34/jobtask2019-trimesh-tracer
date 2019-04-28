/** OptiX SDK 6.0.0 - optix_datatypes.h */
#define RT_DEFAULT_MAX 1.e27f

typedef struct
{
	/** The origin of the ray */
	float3 origin;
	/** The direction of the ray */
	float3 direction;
	/** The ray type associated with this ray */
	unsigned int ray_type;
	/** The min extent associated with this ray */
	float tmin;
	/** The max extent associated with this ray */
	float tmax;
} Ray;

static inline Ray make_Ray(float3 origin, float3 direction, unsigned int ray_type, float tmin, float tmax)
{
	Ray ray;
	
	ray.origin = origin;
	ray.direction = direction;
	ray.ray_type = ray_type;
	ray.tmin = tmin;
	ray.tmax = tmax;
	
	return ray;
}

/** aras-p/jobtask2019-trimesh-tracer */
// general minimum/maximum distances for rays
__constant float kMinT = 0.001f;
__constant float kMaxT = 1.0e7f;
// maximum raytracing recursion depth, i.e. number of light bounces
__constant int kMaxDepth = 10;

// we have one hardcoded directional light, with this direction and color
// we need a compile-time constant
//__constant float3 kLightDir = normalize((float3)(-0.7f,1.0f,0.5f));
__constant float3 kLightDir = (float3)(-0.530668616f, 0.758098066f, 0.379049033f);
__constant float3 kLightColor = (float3)(0.7f,0.6f,0.5f);

static inline bool intersect_triangle_raysegment
(
	const Ray*    ray,
	const float3* p0,
	const float3* p1,
	const float3* p2,
		  float3* n,
		  float*  t,
		  float*  beta,
		  float*  gamma
)
{
	float3 edge0 = (*p1) - (*p0);
	float3 edge1 = (*p2) - (*p1);
	float3 normal = normalize(cross(edge0, edge1));
	float planeOffset = dot(*p0, normal);
	
	float3 _p0 = ray->origin + ray->direction * (kMinT);
	float3 _p1 = ray->origin + ray->direction * (kMaxT);
	
	float offset0 = dot(_p0, normal);
	float offset1 = dot(_p1, normal);
	
	// does the ray segment between tMin & tMax intersect the triangle plane?
	if ((offset0 - planeOffset) * (offset1 - planeOffset) <= 0.0f)
	{
		float _t = kMinT + (kMaxT - kMinT)*(planeOffset - offset0) / (offset1 - offset0);
		float3 p = ray->origin + ray->direction * (_t);
		
		float3 c0 = cross(edge0, p - (*p0));
		float3 c1 = cross(edge1, p - (*p1));
		if (dot(c0, c1) >= 0.f)
		{
			float3 edge2 = (*p0) - (*p2);
			float3 c2 = cross(edge2, p - (*p2));
			if (dot(c1, c2) >= 0.0f)
			{
				*n = normal;
				*t = _t;
				return true;
			}
		}
	}
	
	return false;
}

typedef struct
{
	float4 v0, v1, v2;
} Triangle;

typedef struct
{
	float3 pos;
	float3 normal;
	float t;
} Hit;

typedef struct
{
	float4 origin;
	float4 lowerLeftCorner;
	float4 horizontal;
	float4 vertical;
	float4 u;
	float4 v;
	float4 w;
	
	float lensRadius;
	float pad0, pad1, pad2;
} Camera;

static inline unsigned int XorShift32(unsigned int *state)
{
	unsigned int x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 15;
	*state = x;
	return x;
}

static inline float RandomFloat01(unsigned int *state)
{
	return (XorShift32(state) & 0xFFFFFF) / 16777216.0f;
}

static inline float3 RandomInUnitDisk(unsigned int *state)
{
	float px,py;
	do {
		px = 2.0f*RandomFloat01(state) - 1.0f;
		py = 2.0f*RandomFloat01(state) - 1.0f;
	} while (px*px + py*py >= 1.0f);
	return (float3)(px,py,0.0f);
}

static inline float3 RandomUnitVector(unsigned int *state)
{
	float z = RandomFloat01(state) * 2.0f - 1.0f;
	float a = RandomFloat01(state) * 2.0f * M_PI;
	float r = sqrt(1.0f - z * z);
	float x = r * cos(a);
	float y = r * sin(a);
	return (float3)(x, y, z);
}

static inline Ray GetRay(const Camera* cam, float s, float t, unsigned int *state)
{
	float3 rd = cam->lensRadius * RandomInUnitDisk(state);
	float3 offset = cam->u.xyz * rd.x + cam->v.xyz * rd.y;
	return make_Ray(cam->origin.xyz + offset, normalize(cam->lowerLeftCorner.xyz + s*cam->horizontal.xyz + t*cam->vertical.xyz - cam->origin.xyz - offset), 0, kMinT, kMaxT);
}

// Check all the triangles in the scene for a hit, and return the closest one.
int HitScene(int s_TriangleCount, __global const Triangle* restrict gTriangles, const Ray* r, Hit* outHit)
{
	float hitMinT = kMaxT;
	int hitID = -1;
	for (int i = 0; i < s_TriangleCount; ++i)
	{
		float3 n;
		float  t, beta, gamma;
		
		float3 v0 = gTriangles[i].v0.xyz;
		float3 v1 = gTriangles[i].v1.xyz;
		float3 v2 = gTriangles[i].v2.xyz;
		
		//if (HitTriangle(r, s_Triangles[i], tMin, tMax, hit))
		if (intersect_triangle_raysegment(r, &v0, &v1, &v2, &n, &t, &beta, &gamma))
		{
			if (t < hitMinT)
			{
				hitMinT = t;
				hitID = i;
				
				outHit->t = t;
				outHit->pos = r->origin + r->direction * t;
				outHit->normal = n;
			}
		}
	}
	
	return hitID;
}

// trace a ray into the scene, and return the normal color for it (debug)
static inline float3 TraceNormal(int s_TriangleCount, __global const Triangle* restrict gTriangles, const Ray* r, int* inoutRayCount)
{
	*inoutRayCount = *inoutRayCount + 1;
	
	Hit hit;
	int id = HitScene(s_TriangleCount, gTriangles, r, &hit);
	if (id != -1)
	{
		// ray hits something in the scene: return world normal
		return (float3)(fabs(hit.normal.x), fabs(hit.normal.y), fabs(hit.normal.z));
	}
	else
	{
		// ray does not hit anything: return black
		return (float3)(0.0f, 0.0f, 0.0f);
	}
}

// when a ray "r" has just hit a surface at point "hit", decide what to do about it:
// in our very simple case, we assume the surface is perfectly diffuse, so we'll return:
// - surface albedo ("color") in "attenuation"
// - new random ray for the next light bounce in "scattered"
// - illumination from the directional light in "outLightE"
static inline bool Scatter(int s_TriangleCount, __global const Triangle* restrict gTriangles, const Ray* r, const Hit* hit, float3* attenuation, Ray* scattered, float3* outLightE, unsigned int* rngState, int* inoutRayCount)
{
	*outLightE = (float3)(0,0,0);
	
	// model a perfectly diffuse material:
	
	// random point on unit sphere that is tangent to the hit point
	float3 target = hit->pos + hit->normal + RandomUnitVector(rngState);
	*scattered = make_Ray(hit->pos, normalize(target - hit->pos), 0, kMinT, kMaxT);
	
	// make color slightly based on surface normals
	float3 albedo = hit->normal * 0.0f + (float3)(0.7f,0.7f,0.7f);
	*attenuation = albedo;
	
	// explicit directional light by shooting a shadow ray
	*inoutRayCount = *inoutRayCount + 1;
	
	Hit lightHit;
	Ray lightRay = make_Ray(hit->pos, kLightDir, 0, kMinT, kMaxT);
	int id = HitScene(s_TriangleCount, gTriangles, &lightRay, &lightHit);
	if (id == -1)
	{
		// ray towards the light did not hit anything in the scene, so
		// that means we are not in shadow: compute illumination from it
		float3 rdir = r->direction;
		float3 nl = dot(hit->normal, rdir) < 0 ? hit->normal : -hit->normal;
		*outLightE += albedo * kLightColor * (fmax(0.0f, dot(kLightDir, nl)));
	}
	
	return true;
}

// trace a ray into the scene, and return the final color for it
static float3 TraceIterative(int s_TriangleCount, __global const Triangle* restrict gTriangles, const Ray* r, unsigned int* rngState, int* inoutRayCount)
{
	Ray currRay;
	currRay.origin = r->origin;
	currRay.direction = r->direction;
	currRay.ray_type = r->ray_type;
	currRay.tmin = r->tmin;
	currRay.tmax = r->tmax;
	
	float3 currAttenuation = (float3)(1.0f, 1.0f, 1.0f);
	float3 currLightE = (float3)(0.0f, 0.0f, 0.0f);
	
	for (int currDepth = 0; currDepth < kMaxDepth; ++currDepth)
	{
		*inoutRayCount = *inoutRayCount + 1;
		
		Hit hit;
		int id = HitScene(s_TriangleCount, gTriangles, &currRay, &hit);
		if (id != -1)
		{
			// ray hits something in the scene
			Ray scattered;
			float3 attenuation;
			float3 lightE;
			if (Scatter(s_TriangleCount, gTriangles, &currRay, &hit, &attenuation, &scattered, &lightE, rngState, inoutRayCount))
			{
				currLightE += lightE * currAttenuation;
				currAttenuation *= attenuation;
				
				currRay.origin = scattered.origin;
				currRay.direction = scattered.direction;
				currRay.ray_type = scattered.ray_type;
				currRay.tmin = scattered.tmin;
				currRay.tmax = scattered.tmax;
			}
			else
			{
				// reached recursion limit, or surface fully absorbed the ray: return black
				return (float3)(0.0f, 0.0f, 0.0f);
			}
		}
		else
		{
			// ray does not hit anything: return illumination from the sky (just a simple gradient really)
			float3 unitDir = currRay.direction;
			float t = 0.5f*(unitDir.y + 1.0f);
			float3 c = ((1.0f - t)*(float3)(1.0f, 1.0f, 1.0f) + t * (float3)(0.5f, 0.7f, 1.0f)) * 0.5f;
			return currLightE + currAttenuation * c;
		}
	}
	
	return (float3)(0.0f, 0.0f, 0.0f); // exceeded recursion
}

__kernel void trace
(
	int samplesPerPixel,
	int width, int height, float invWidth, float invHeight,
	int triangleNum,
	__global const Triangle* restrict gTriangles,
	__global unsigned int* restrict gRngState,
	__global const Camera* restrict gCamera,
	__global uchar4* restrict image
)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	
	int rayCount = 0;
	
	unsigned int rngState = gRngState[x + y * width];
	
	Camera cam = gCamera[0];
	
	float3 col = (float3)(0.0f, 0.0f, 0.0f);
	// we'll trace N slightly jittered rays for each pixel, to get anti-aliasing, loop over them here
	for(int s=0; s < samplesPerPixel; ++s)
	{
		// get a ray from camera, and trace it
		float u = (float)(x + RandomFloat01(&rngState)) * invWidth;
		float v = (float)(y + RandomFloat01(&rngState)) * invHeight;
		Ray r = GetRay(&cam, u, v, &rngState);
		
		//col += TraceNormal(triangleNum, gTriangles, &r, &rayCount);
		col += TraceIterative(triangleNum, gTriangles, &r, &rngState, &rayCount);
	}
	col *= 1.0f / (float)(samplesPerPixel);
	
	// simplistic "gamma correction" by just taking a square root of the final color
	col = sqrt(col);
	
	// our image is bytes in 0-255 range, turn our floats into them here and write into the image
	col = clamp(col, (float3)(0.0f,0.0f,0.0f), (float3)(1.0f,1.0f,1.0f));
	image[x + y * width] = (uchar4)(col.x * 255, col.y * 255, col.z * 255, 255);
}