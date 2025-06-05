/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
 */

 #include "./common.glsl"
 #iChannel0 "self"

 #define SCENE 1

bool hit_world(Ray r, float tmin, float tmax, inout HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;

    #if SCENE == 0       //Shirley Weekend scene

        if(hit_quad(createQuad(vec3(-10.0, -0.05, 10.0), vec3(10.0, -0.05, 10.0), vec3(10.0, -0.05, -10.0), vec3(-10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2));
        }

        if(hit_sphere(createSphere(vec3(-4.0, 1.0, 0.0), 1.0), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
            rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
        }

        if(hit_sphere(createSphere(vec3(4.0, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), -0.5),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0, 0.9, 0.9), 1.5, 0.0);
        }
            
        int numxy = 5;
        
        for(int x = -numxy; x < numxy; ++x)
        {
            for(int y = -numxy; y < numxy; ++y)
            {
                float fx = float(x);
                float fy = float(y);
                float seed = fx + fy / 1000.0;
                vec3 rand1 = hash3(seed);
                vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
                float chooseMaterial = rand1.z;
                if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
                {
                    if(chooseMaterial < 0.3)
                    {
                        vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                        // diffuse
                        if(hit_movingSphere(createMovingSphere(center, center1, 0.2, 0.0, 1.0),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.5)
                    {
                        // diffuse
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.7)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                        }
                    }
                    else if(chooseMaterial < 0.9)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                        }
                    }
                    else
                    {
                        // glass (Dielectric)
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDielectricMaterial(hash3(seed), 1.33, 0.0);
                        }
                    }
                }
            }
        }
    #elif SCENE == 1 //from https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/

        // diffuse floor
        
            vec3 A = vec3(-25.0f, -12.5f, 10.0f);
            vec3 B = vec3( 25.0f, -12.5f, 10.0f);
            vec3 C = vec3( 25.0f, -12.5f, -5.0f);
            vec3 D = vec3(-25.0f, -12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }

        //stripped background
        {
            vec3 A = vec3(-25.0f, -10.5f, -5.0f);
            vec3 B = vec3( 25.0f, -10.5f, -5.0f);
            vec3 C = vec3( 25.0f, -1.5f, -5.0f);
            vec3 D = vec3(-25.0f, -1.5f, -5.0f);
        
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                float shade = floor(mod(rec.pos.x, 1.0f) * 2.0f);
                rec.material = createDiffuseMaterial(vec3(shade));
            }
        }

        // ceiling piece above light
        
        {
            vec3 A = vec3(-7.5f, 12.5f, 5.0f);
            vec3 B = vec3( 7.5f, 12.5f, 5.0f);
            vec3 C = vec3( 7.5f, 12.5f, -5.0f);
            vec3 D = vec3(-7.5f, 12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }
        }    
       
        // light
        
        {
            vec3 A = vec3(-5.0f, 12.3f,  2.5f);
            vec3 B = vec3( 5.0f, 12.3f,  2.5f);
            vec3 C = vec3( 5.0f, 12.3f,  -2.5f);
            vec3 D = vec3(-5.0f, 12.3f,  -2.5f);

             if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.9f) * 20.0f;
            }
        }
 
        const int c_numSpheres = 7;
        for (int sphereIndex = 0; sphereIndex < c_numSpheres; ++sphereIndex)
        {
            vec3 center = vec3(-18.0 + 6.0 * float(sphereIndex), -8.0, 0.0);
            if(hit_sphere(createSphere(center, 2.8),r,tmin,rec.t,rec))
            {
                hit = true;
                float r = float(sphereIndex) / float(c_numSpheres-1) * 0.1f;
                rec.material = createDielectricMaterial(vec3(0.0, 0.5, 1.0), 1.1, r);
            }
        }

    #elif SCENE == 2
        // Scene with just a sphere. 
        // This is a simple test scene to verify basic ray-sphere intersection.
        if(hit_sphere(createSphere(vec3(0.0, 0.0, -5.0), 1.0), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.8, 0.2, 0.2));
        }
    #elif SCENE == 3
    #endif

    return hit;
}

vec3 directlighting(pointLight pl, Ray r, HitRecord rec) {
    vec3 colorOut = vec3(0.0);
    
    // Vector from hit point to light
    vec3 lightDir = normalize(pl.pos - rec.pos);
    
    // Distance to light (for attenuation)
    float lightDistance = length(pl.pos - rec.pos);
    
    // Check if light is on the same side as the surface normal
    float NdotL = dot(rec.normal, lightDir);
    if (NdotL <= 0.0) {
        return vec3(0.0); // Light is behind the surface
    }
    
    // Shadow ray - check if there's an object between hit point and light
    Ray shadowRay;
    shadowRay.o = rec.pos + rec.normal * 0.001; // Offset to avoid self-intersection
    shadowRay.d = lightDir;
    
    HitRecord shadowRec;
    // If shadow ray hits something before reaching the light, we're in shadow
    if (hit_world(shadowRay, 0.001, lightDistance - 0.001, shadowRec)) {
        return vec3(0.0); // In shadow
    }
    
    // Light attenuation (inverse square law)
    float attenuation = 1.0 / (1.0 + 0.1 * lightDistance + 0.01 * lightDistance * lightDistance);
    
    // Calculate lighting based on material type
    if (rec.material.type == MT_DIFFUSE) {
        // Lambertian diffuse
        vec3 diffuse = rec.material.albedo * pl.color * NdotL * attenuation;
        colorOut += diffuse;
    }
    else if (rec.material.type == MT_METAL) {
        // Blinn-Phong for metals
        vec3 viewDir = normalize(-r.d);
        vec3 halfDir = normalize(lightDir + viewDir);
        float NdotH = max(dot(rec.normal, halfDir), 0.0);
        
        // Diffuse component
        vec3 diffuse = rec.material.albedo * pl.color * NdotL;
        
        // Specular component
        float shininess = mix(32.0, 128.0, 1.0 - rec.material.roughness);
        vec3 specular = rec.material.albedo * pl.color * pow(NdotH, shininess);
        
        colorOut += (diffuse + specular) * attenuation;
    }
    else if (rec.material.type == MT_DIELECTRIC) {
        // For glass/dielectric materials, mainly handle transmission
        // Simple approach: some diffuse reflection
        vec3 diffuse = rec.material.albedo * pl.color * NdotL * 0.1; // Reduced contribution
        colorOut += diffuse * attenuation;
    }
    
    return colorOut;
}

#define MAX_BOUNCES 10

vec3 rayColor(Ray r)
{
    HitRecord rec;
    vec3 col = vec3(0,0,0);
    vec3 throughput = vec3(1.0);
    
    for(int i = 0; i < MAX_BOUNCES; ++i)
    {
        if(hit_world(r, 0.001, 10000.0, rec))
        {
            // Add emissive contribution (for light sources)
            col += rec.material.emissive * throughput;
            
            // Calculate direct lighting with 3 white point lights
            col += directlighting(createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
            col += directlighting(createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
            col += directlighting(createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
            col = vec3(1,1,1);

            // Calculate secondary ray and update throughput
            Ray scatterRay;
            vec3 atten;

            // if(scatter(r, rec, atten, scatterRay))
            // {
            //     // Update throughput with material attenuation
            //     throughput *= atten;
                
            //     // Continue with the scattered ray
            //     r = scatterRay;
                
            //     // Russian Roulette for performance optimization (optional)
            //     // Terminate rays that contribute very little to the final image
            //     if(i > 3) // Start Russian Roulette after a few bounces
            //     {
            //         float maxComponent = max(throughput.r, max(throughput.g, throughput.b));
            //         if(maxComponent < 0.1) // Threshold for termination
            //         {
            //             float continueProbability = maxComponent;
            //             if(hash1(gSeed) > continueProbability)
            //             {
            //                 break; // Terminate ray
            //             }
            //             // Boost throughput to maintain energy conservation
            //             throughput /= continueProbability;
            //         }
            //     }
            // }
            // else
            // {
            //     // Material absorbed the ray (no scattering occurred)
            //     break;
            // }
        }
        else // Ray missed all objects - hit background/environment
        {
            // Simple sky gradient background
            
            float t = 0.5 * (r.d.y + 1.0); // Normalize y to [0,1]
            vec3 skyColor = mix(vec3(1.0, 1.0, 1.0), vec3(0.5, 0.7, 1.0), t);
            col += throughput * skyColor;  
            col = vec3(t-0.5,0,0);

            break;
        }
    }
    
    return col;
    // return vec3(0.5 * (r.d + 1.0)); // Should show a colorful gradient
}

#define MAX_SAMPLES 10000.0

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    vec2 mouse = iMouse.xy / iResolution.xy;
    mouse.x = mouse.x * 2.0 - 1.0;

    vec3 camPos = vec3(mouse.x * 10.0, mouse.y * 5.0, 8.0);
    vec3 camTarget = vec3(0.0, 0.0, -1.0);
    float fovy = 60.0;
    float aperture = 0.0;
    float distToFocus = 1.0;
    float time0 = 0.0;
    float time1 = 1.0;
    Camera cam = createCamera(
        camPos,
        camTarget,
        vec3(0.0, 1.0, 0.0),    // world up vector
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

    //usa-se o 4 canal de cor para guardar o numero de samples e não o iFrame pois quando se mexe o rato faz-se reset

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));
    // vec3 color = vec3(1,1,0);

    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.0;
    color = mix(prevLinear, color, 1.0/w);
    gl_FragColor = vec4(toGamma(color), w);
}
