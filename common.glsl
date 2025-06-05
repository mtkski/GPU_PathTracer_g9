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

float hash1(inout float seed) 
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) 
{
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

Ray getRay(Camera cam, vec2 pixel_sample) 
{
    // Lens sample 
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    // Eye offset 
    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;
    
    // OPTIMISATION: Calc directly on vp
    // pixel_sample is in [0,1]x[0,1], we center directly
    vec2 centered_pixel = pixel_sample - 0.5; // [-0.5, 0.5]
    
    // Point on the vp in world space coords
    vec3 viewport_center = cam.eye - cam.n * cam.planeDist;
    vec3 viewport_point = viewport_center + 
                         centered_pixel.x * cam.width * cam.u + 
                         centered_pixel.y * cam.height * cam.v;
    
    vec3 ray_direction = normalize(viewport_point - eye_offset);
    
    return createRay(eye_offset, ray_direction, time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIELECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for Dielectric
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

Material createDielectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIELECTRIC;
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
    float eta_i = 1.0; //ior_1 = 1.0;
    // refIdx = eta_t
    float r0 = (eta_i - refIdx) / (eta_i + refIdx);
    r0 = r0 * r0;
    float Kr = r0 + (eta_i - r0) * pow(eta_i - cosine, 5.0);
    return Kr;
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        //scatter direction
        vec3 scatterDir = rec.normal + randomUnitVector(gSeed);
        scatterDir = normalize(scatterDir); //ensure direction is normalized
        //if scatterDir is too close to zero, use rec.normal
        if(dot(scatterDir, rec.normal) < 0.0) scatterDir = rec.normal; //ensure we don't scatter into the back side
        rScattered = createRay(rec.pos, scatterDir, rIn.t);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;     
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections
        atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIELECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt;
        float cosine;

        if(dot(rIn.d, rec.normal) > 0.0) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            //cosine = refraction cosine for schlick;
            cosine = dot(rIn.d, outwardNormal);
            //atten = apply Beers law by using rec.material.refractColor;

        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rIn.d, rec.normal); 
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray

        float reflectProb;

        //if no total reflection  reflectProb = schlick(cosine, rec.material.refIdx);  
        //else reflectProb = 1.0;

        if( hash1(gSeed) < reflectProb)  //Reflection
        // rScattered = calculate reflected ray
         
        //else  //Refraction
        // rScattered = calculate refracted ray
          
        

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
    vec3 edgeBA = t.b - t.a;
    vec3 edgeCA = t.c - t.a;
    vec3 P = cross(r.d, edgeCA);
    vec3 edgeCB = t.c - t.b;
    vec3 normal = normalize(cross(edgeBA, edgeCB));

    float det = dot(edgeBA, P);

    if(det < epsilon && det > -epsilon) return false; //ray is parallel to triangle

    float invDet = 1.0 / det;
    vec3 T = r.o - t.a;
    vec3 Q = cross(T, edgeBA);
    
    float tHit = dot(edgeCA, Q) * invDet;
    // float tHit = 0.2;

    if(tHit < tmin  ) return false; //hit outside of range

    float beta = dot(T, P) * invDet;

    if(beta < 0.0 || beta > 1.0) return false;

    float gamma = dot(r.d, Q) * invDet;

    if(gamma < 0.0 || beta + gamma > 1.0) return false;

    if(tHit < tmax && tHit > tmin) //hit inside range
    {
        rec.t = tHit;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        return true;
    }

    return false;  //hit outside of range
}


struct Quad {vec3 a; vec3 b; vec3 c; vec3 d; };

Quad createQuad(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{
    Quad q;
    q.a = v0; q.b = v1; q.c = v2; q.d = v3;
    return q;
}

bool hit_quad(Quad q, Ray r, float tmin, float tmax, out HitRecord rec)
{
    if(hit_triangle(createTriangle(q.a, q.b, q.c), r, tmin, rec.t, rec)) return true;
    else if(hit_triangle(createTriangle(q.a, q.c, q.d), r, tmin, rec.t, rec)) return true;
    else return false;  
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
	//Program it
    vec3 moving_center = mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * (time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0);
    return moving_center;
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    //calculate a valid t and normal
    vec3 OC = s.center - r.o;
    float b = dot(OC, r.d);
    float c = dot(OC, OC) - s.radius * s.radius;
    if(c > 0.0f){
        if(b < 0.0f) return false; //ray is outside sphere and pointing away
    }
    float discriminant = b * b - c;
    if(discriminant < 0.0f) return false; //no intersection
    float sqrtDiscriminant = sqrt(discriminant);
    if(c > 0.0f) {
        // ray is outside sphere
        float t = b - sqrtDiscriminant; // first intersection
        if(t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            vec3 HitPoint = rec.pos;
            rec.normal = normalize(HitPoint - s.center); // normal is from center to hit point
            return true;
    }
    } else {
        // ray is inside sphere
        float t = b + sqrtDiscriminant; // first intersection
        if(t < tmax && t > tmin) {
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
            vec3 HitPoint = rec.pos;
            rec.normal = normalize(HitPoint - s.center); // normal is from center to hit point
            return true;
        }
    }
    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    float B, C, delta;
    bool outside;
    float t;


     //INSERT YOUR CODE HERE
     //Calculate the moving center
    vec3 center = center(s, r.t);
    vec3 OC = center - r.o;
    B = dot(OC, r.d);
    C = dot(OC, OC) - s.radius * s.radius;
    delta = B * B - C;
    if(delta < 0.0f) return false; //no intersection
    float sqrtDelta = sqrt(delta);

    //calculate a valid t and normal
    if(C > 0.0f) {
        // ray is outside sphere
        t = B - sqrtDelta; // first intersection
        outside = true;
    } else {
        // ray is inside sphere
        t = B + sqrtDelta; // first intersection
        outside = false;
    }
    vec3 normal = normalize(r.o + r.d * t - center); // normal is from center to hit point
    if(!outside) {
        normal = -normal; // if inside, normal points in the opposite direction
    }
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        return true;
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