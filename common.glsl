/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction
    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;

    float x = cam.width * (pixel_sample.x / iResolution.x - 0.5) * cam.focusDist;
    float y = cam.height * (pixel_sample.y / iResolution.y - 0.5) * cam.focusDist;
    float z = cam.planeDist * cam.focusDist;

    vec3 p = vec3(x, y, -z) - vec3(ls, 0.0);

    vec3 ray_direction = cam.u * p.x + cam.v * p.y + cam.n * p.z;

    return createRay(eye_offset, normalize(ray_direction), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float refIdx)
{
    //INSERT YOUR CODE HERE
    float r0 = pow((refIdx - 1.0) / (refIdx + 1.0), 2.0);
    float x = 1.0 - cosine;
    float ret = r0 + (1.0 - r0) * pow(x, 5.0);

    return ret;
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{

    if(rec.material.type == MT_DIFFUSE)
    {
        //INSERT CODE HERE,
        vec3 rDir = rec.normal + normalize(randomUnitVector(gSeed));
        rScattered = createRay(rec.pos, rDir, rec.t);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections
        vec3 fuzzyReflcetionsDir = reflect(rIn.d, rec.normal) + rec.material.roughness * randomInUnitSphere(gSeed);
        rScattered = createRay(rec.pos + rec.normal * epsilon, normalize(fuzzyReflcetionsDir), rIn.t);
        atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt, cosine;

        bool inside = dot(rIn.d, rec.normal) > 0.0;

        if(inside) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            cosine = dot(rIn.d, rec.normal) * rec.material.refIdx;
            atten = exp(-rec.material.refractColor * rec.t);
        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rIn.d, rec.normal); 
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray
        float reflectProb;
        vec3 v = -rIn.d;
        vec3 vt = (outwardNormal * dot(v, outwardNormal)) - v;
        float sinI = length(vt);
        float sinT = niOverNt * sinI;

        float total_reflection = sinT * sinT;

        if (total_reflection < 1.0) {
            reflectProb = schlick(cosine, rec.material.refIdx);
        } else {
            reflectProb = 1.0;
        }

        if (hash1(gSeed) < reflectProb) {
            //Reflection
            vec3 rOrigin = inside? rec.pos - rec.normal * epsilon : rec.pos + rec.normal * epsilon;
            rScattered = createRay(rOrigin, reflect(rIn.d, rec.normal), rIn.t);
            // atten *= vec3(reflectProb); not necessary since we are only scattering reflectProb rays and not all reflected rays
        } else {
            //Refraction
            vec3 rOrigin = inside? rec.pos + rec.normal * epsilon : rec.pos - rec.normal * epsilon;
            rScattered = createRay(rOrigin, refract(rIn.d, rec.normal, niOverNt), rIn.t);
            // atten *= vec3(1.0 - reflectProb); not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
        }

        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    vec3 AB = t.b - t.a;
    vec3 AC = t.c - t.a;
    vec3 rOA = r.o - t.a;

    vec3 normal = normalize(cross(AB, AC));
    float temp = dot(-r.o, normal) / dot(normal, r.d);
    vec3 P = r.o + r.d * temp;

    float s1 = t.c.y - t.a.x;
    float s2 = t.c.x - t.a.x;
    float s3 = t.b.y - t.a.y;
    float s4 = P.y - t.a.y;

    float w1 = (t.a.x * s1 + s4 * s2 - P.x * s1) / (s3 * s2 - (t.b.x-t.a.x) * s1);
    float w2 = (s4 - w1 * s3) / s1;

    if (w1 < 0.0 && w2 < 0.0 && w1 + w2 > 1.0) return false;

    //calculate a valid t and normal
    if(temp < tmax && temp > tmin)
    {
        rec.t = temp;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }

    return false;
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    return mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * time;
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    float sR2 = s.radius * s.radius;
    vec3 OC = s.center - r.o;
    float A = dot(r.d, r.d);
    float B = dot(OC, r.d);
    float C2 = dot(OC, OC) - sR2;
    float discriminant = B * B - A * C2;

    if (discriminant > 0.0) {
        float tH = sqrt(discriminant);
        float t = (B - tH) / A;
        
        if (t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = normalize(rec.pos - s.center);
            //if (s.radius < 0.0) rec.normal = -rec.normal; //inside sphere
            return true;
        }

        t = (B + tH) / A;
        if (t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = normalize(rec.pos - s.center);
            //if (s.radius < 0.0) rec.normal = -rec.normal; //inside sphere
            return true;
        }
    }

    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    vec3 center = center(s, r.t);
    float sR2 = s.radius * s.radius;
    vec3 OC = center - r.o;
    float A = dot(r.d, r.d);
    float B = dot(OC, r.d);
    float C2 = dot(OC, OC) - sR2;
    float discriminant = B * B - A * C2;

    if (discriminant > 0.0) {
        float tH = sqrt(discriminant);
        float t = (B - tH) / A;
        
        if (t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = normalize(rec.pos - center);
            return true;
        }

        t = (B + tH) / A;
        if (t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = normalize(rec.pos - center);
            return true;
        }
    }

    return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}