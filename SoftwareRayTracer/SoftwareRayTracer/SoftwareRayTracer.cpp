#define _USE_MATH_DEFINES 
#include <iostream>
#include <map>
#include <vector>

//#include <SDL.h>
#include <Windows.h>
#undef main

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "bmp.h"
#include "obj.h"
#include "verts.h"
#include "BVH.h"

#define PIXEL_W 1200
#define PIXEL_H 800

//SDL_Event event;
//SDL_Window* window;
//SDL_Renderer* renderer;

/*
float verts[] =
{
    //pos               col                 
    -.5f, -.5f, 2.f,    1.f, 0.f, 0.f,      //tl
    .5f, -.5f, 2.f,     1.f, 0.f, 0.f,      //tr
    .5f, .5f, 2.f,    1.f, 0.f, 0.f,      //br
    -.5f, -.5f, 2.f,    0.f, 1.f, 0.f,      //tl
    .5f, .5f, 2.f,     0.f, 1.f, 0.f,      //tr
    -.5f, .5f, 2.f,    0.f, 1.f, 0.f,      //br
};
*/

struct payload
{
    vec3 colour;
    int  depth;
    int  rayHitSky;
    vec3 rayOrigin;
    vec3 rayDir;
    uint rngState;
    vec3 hitNormal;
};


glm::vec3 DoNothing(triangle* tri, int depth, glm::vec3 p, glm::vec3 dir);
glm::vec3 CalculateColourFlat(triangle* tri, int depth, glm::vec3 p, glm::vec3 dir, payload* prd);

vec3 skyColour(vec3 rayDir)
{
    vec3 colour;
    // +y in world space is up, so:
    if (rayDir.y > 0.0f)
    {
        colour = mix(vec3(1.0f), vec3(0.25f, 0.5f, 1.0f), rayDir.y);
    }
    else
    {
        colour = vec3(0.03f);
    }
    return colour;
}

glm::vec3 bkgd = glm::vec3(1.f, 1.f, 1.f);

const float atten = 1.f;

int max_recursion_depth = 1;

const int use_bvh = 1;

int samples_per_pixel = 16;
int tracedSegments = 8;

const int bvh_width = 2;

vector<Object> objs;
std::vector<triangle> tris;
BVH_node g_BVH;
float pixelBuffer[PIXEL_W * PIXEL_H * 3];

glm::vec3 eye = glm::vec3(0.f, .05f, 2.15f);
float l = -1.f;
float r = 1.f;
float t = 1.f;
float b = -1.f;
float n = 1.f;
float f = 10.f;

glm::vec3 light_pos(1.f, 2.f, 1.f);


void AppendTriangles(std::vector<triangle>* io, vector<Object> in_objs)
{
    for (auto& obj : in_objs)
    {
        io->insert(io->end(), obj.tris.begin(), obj.tris.end());
    }
}

std::vector<triangle> AssemblePrimitives(float* verts, int n_verts)
{
    int n_tris = n_verts / 3;
    std::vector<triangle> ret;

    for (int tc = 0; tc < n_tris; tc++)
    {
        triangle t;

        t.v1.pos = glm::vec3(verts[(tc * 30) + 0], verts[(tc * 30) + 1], verts[(tc * 30) + 2]);
//        t.v1.pos.x = verts[(tc * 30) + 0];
//        t.v1.pos.y = verts[(tc * 30) + 1];
//        t.v1.pos.z = verts[(tc * 30) + 2];
        t.v1.col = glm::vec3(verts[(tc * 30) + 3], verts[(tc * 30) + 4], verts[(tc * 30) + 5]);
//        t.v1.col.x = verts[(tc * 30) + 3];
//        t.v1.col.y = verts[(tc * 30) + 4];
//        t.v1.col.z = verts[(tc * 30) + 5];
        t.v1.nor = glm::vec3(verts[(tc * 30) + 7], verts[(tc * 30) + 8], verts[(tc * 30) + 9]);
//        t.v1.nor.x = verts[(tc * 30) + 7];
//        t.v1.nor.y = verts[(tc * 30) + 8];
//        t.v1.nor.z = verts[(tc * 30) + 9];

        t.v2.pos = glm::vec3(verts[(tc * 30) + 10], verts[(tc * 30) + 11], verts[(tc * 30) + 12]);
//        t.v2.pos.x = verts[(tc * 30) + 10];
//        t.v2.pos.y = verts[(tc * 30) + 11];
//        t.v2.pos.z = verts[(tc * 30) + 12];
        t.v2.col = glm::vec3(verts[(tc * 30) + 13], verts[(tc * 30) + 14], verts[(tc * 30) + 15]);
//        t.v2.col.x = verts[(tc * 30) + 13];
//        t.v2.col.y = verts[(tc * 30) + 14];
//        t.v2.col.z = verts[(tc * 30) + 15];
        t.v2.nor = glm::vec3(verts[(tc * 30) + 17], verts[(tc * 30) + 18], verts[(tc * 30) + 19]);
//        t.v2.nor.x = verts[(tc * 30) + 17];
//        t.v2.nor.y = verts[(tc * 30) + 18];
//        t.v2.nor.z = verts[(tc * 30) + 19];

        t.v3.pos = glm::vec3(verts[(tc * 30) + 20], verts[(tc * 30) + 21], verts[(tc * 30) + 22]);
//        t.v3.pos.x = verts[(tc * 30) + 20];
//        t.v3.pos.y = verts[(tc * 30) + 21];
 //       t.v3.pos.z = verts[(tc * 30) + 22];
        t.v3.col = glm::vec3(verts[(tc * 30) + 23], verts[(tc * 30) + 24], verts[(tc * 30) + 25]);
//        t.v3.col.x = verts[(tc * 30) + 23];
//        t.v3.col.y = verts[(tc * 30) + 24];
 //       t.v3.col.z = verts[(tc * 30) + 25];
        t.v3.nor = glm::vec3(verts[(tc * 30) + 27], verts[(tc * 30) + 28], verts[(tc * 30) + 29]);
//        t.v3.nor.x = verts[(tc * 30) + 27];
//        t.v3.nor.y = verts[(tc * 30) + 28];
 //       t.v3.nor.z = verts[(tc * 30) + 29];

        t.reflect = verts[(tc * 30) + 26] ? true : false;
        t.primitiveID = tc;
        ret.push_back(t);
    }

    return ret;
}

glm::vec3 GetPixelInViewSpace(int pixel_x, int pixel_y, int W, int H, float l, float r, float t, float b, float n, float f)
{
    glm::vec3 pvs;

    float xr = float(pixel_x) / float(W);
    float yr = float(pixel_y) / float(H);

    pvs.x = l + ((r - l) * xr);
    pvs.y = b + ((t - b) * yr);
    pvs.z = n;

    return pvs;
}

glm::vec3 GetRandomPixelPosInViewSpace(int pixel_x, int pixel_y, int W, int H, float l, float r, float t, float b, float n, float f)
{
    glm::vec3 pvs;

    float xr = float(pixel_x) / float(W);
    float yr = float(pixel_y) / float(H);

    float randx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    randx *= 2.f;
    randx -= 1.f;
    randx *= 1.f/float(W);
    float randy = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    randy *= 2.f;
    randy -= 1.f;
    randy *= 1.f / float(H);


    pvs.x = l + ((r - l) * xr);
    pvs.x += randx;
    pvs.y = b + ((t - b) * yr);
    pvs.y += randy;
    pvs.z = n;

    return pvs;
}

/*
bool ComputeBarycentricCoordinates(float x_ndc, float y_ndc, triangle& t, float& alpha, float& beta, float& gamma)
{
    float bc_edge = (t.v2.pos.y - t.v1.pos.y) * x_ndc + (t.v1.pos.x - t.v2.pos.x) * y_ndc + (t.v2.pos.x * t.v1.pos.y) - (t.v1.pos.x * t.v2.pos.y);
    float bc_point = (t.v2.pos.y - t.v1.pos.y) * t.v3.pos.x + (t.v1.pos.x - t.v2.pos.x) * t.v3.pos.y + (t.v2.pos.x * t.v1.pos.y) - (t.v1.pos.x * t.v2.pos.y);
    alpha = bc_edge / bc_point;

    float ac_edge = (t.v3.pos.y - t.v1.pos.y) * x_ndc + (t.v1.pos.x - t.v3.pos.x) * y_ndc + (t.v3.pos.x * t.v1.pos.y) - (t.v1.pos.x * t.v3.pos.y);
    float ac_point = (t.v3.pos.y - t.v1.pos.y) * t.v2.pos.x + (t.v1.pos.x - t.v3.pos.x) * t.v2.pos.y + (t.v3.pos.x * t.v1.pos.y) - (t.v1.pos.x * t.v3.pos.y);
    beta = ac_edge / ac_point;

    float ab_edge = (t.v3.pos.y - t.v2.pos.y) * x_ndc + (t.v2.pos.x - t.v3.pos.x) * y_ndc + (t.v3.pos.x * t.v2.pos.y) - (t.v2.pos.x * t.v3.pos.y);
    float ab_point = (t.v3.pos.y - t.v2.pos.y) * t.v1.pos.x + (t.v2.pos.x - t.v3.pos.x) * t.v1.pos.y + (t.v3.pos.x * t.v2.pos.y) - (t.v2.pos.x * t.v3.pos.y);
    gamma = ab_edge / ab_point;

    if (alpha > 0.f && alpha < 1.f &&
        beta > 0.f && beta < 1.f &&
        gamma > 0.f && gamma < 1.f)
    {
        return true;
    }

    return false;
}
*/
/*
float sign(vec3 p1, vec3 p2, vec3 p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool PointInTriangle(vec3 pt, vec3 v1, vec3 v2, vec3 v3)
{
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(pt, v1, v2);
    d2 = sign(pt, v2, v3);
    d3 = sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}
*/
/*
bool PointInTriangle(vec3 pt, vec3 v1, vec3 v2, vec3 v3)
{
    float edge1 = (v2.y - v1.y) * pt.x + (v1.x - v2.x) * pt.y + (v2.x * v1.y) - (v1.x * v2.y);
    float edge2 = (v3.y - v2.y) * pt.x + (v2.x - v3.x) * pt.y + (v3.x * v2.y) - (v2.x * v3.y);
    float edge3 = (v1.y - v3.y) * pt.x + (v3.x - v1.x) * pt.y + (v1.x * v3.y) - (v3.x * v1.y);

    if (edge1 < 0.f && edge2 < 0.f && edge3 < 0.f)
        return true;

    return false;
}
*/

bool PointInTriangle(glm::vec3 pt, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
{
    glm::vec3 V12 = v2 - v1;
    glm::vec3 V13 = v3 - v1;
    glm::vec3 V1p = pt - v1;
    glm::vec3 C11 = cross(V12, V1p);
    glm::vec3 C12 = cross(V12, V13);
    float d1 = dot(C11, C12);

    glm::vec3 V23 = v3 - v2;
    glm::vec3 V21 = v1 - v2;
    glm::vec3 V2p = pt - v2;
    glm::vec3 C21 = cross(V23, V2p);
    glm::vec3 C22 = cross(V23, V21);
    float d2 = dot(C21, C22);

    glm::vec3 V31 = v1 - v3;
    glm::vec3 V32 = v2 - v3;
    glm::vec3 V3p = pt - v3;
    glm::vec3 C31 = cross(V31, V3p);
    glm::vec3 C32 = cross(V31, V32);
    float d3 = dot(C31, C22);

    if (d1 > 0 && d2 > 0 && d3 > 0)
        return true;

    return false;
}


float RayTriangleIntersection(glm::vec3 o, glm::vec3 dir, triangle *tri, glm::vec3 &point)
{
    //on triangle plane?
    glm::vec3 V1 = tri->v2.pos - tri->v1.pos;
    glm::vec3 V2 = tri->v3.pos - tri->v1.pos;
    glm::vec3 n = cross(V2, V1);
    glm::vec3 N = normalize(n);

    glm::vec3 p = tri->v1.pos - o;
    float d1 = dot(p, N);
    float d2 = dot(dir, N);
    if (d2 == 0.f)
    {
        return FLT_MAX;
    }
    float t1 = d1 / d2;

    if (t1 < 0.001f)
        return FLT_MAX;

    glm::vec3 V3 = dir * t1;
    glm::vec3 pop = o + V3;


    //inside triangle?
    float alpha, beta, gamma;
    if (PointInTriangle(pop, tri->v1.pos, tri->v2.pos, tri->v3.pos))
    {
        point = pop;
        return t1;
    }

    return FLT_MAX;
}



//void RayTrianglesIntersection(glm::vec3 o, glm::vec3 dir, float& t, glm::vec3& io_col, int depth, triangle* from_tri, closest_hit p_hit);
//void RayBVHIntersection(glm::vec3 o, glm::vec3 dir, BVH_node* inBVH, float& t, glm::vec3& ioCol, int depth, triangle* from_tri, closest_hit p_hit, payload *prd);




glm::vec3 Diffuse(glm::vec3 col, glm::vec3 lightDir, glm::vec3 normal)
{
    float dotNL = std::max(dot(normal, lightDir), 0.f);
    glm::vec3 ret = col * dotNL;
    return ret;
}

glm::vec3 DoNothing(triangle* tri, int depth, glm::vec3 p, glm::vec3 dir)
{
    return glm::vec3(0.f);
}
/*
glm::vec3 CalculateColourWhitted(triangle *tri, int depth, glm::vec3 p, glm::vec3 dir)
{
    glm::vec3 refl_col(0.f);
    float t = -1;

    glm::vec3 amb(0.f), diff(0.f);

    amb = .1f * tri->v1.col;

    
    //shadow
    glm::vec3 dummy;
    glm::vec3 lightDir = glm::normalize(light_pos - p);
    if (use_bvh)
        RayBVHIntersection(p, lightDir, &g_BVH, t, dummy, 0, tri, DoNothing);
    else
        RayTrianglesIntersection(p, lightDir, t, dummy, 0, tri, DoNothing);
    if (t < 0)
    {
        diff = Diffuse(tri->v1.col, lightDir, tri->v1.nor);
    }
    

    //reflection
    if (tri->reflect && depth < max_recursion_depth)
    {
        glm::vec3 refl = glm::reflect(dir, tri->v1.nor);
        if(use_bvh)
           RayBVHIntersection(p, refl, &g_BVH, t, refl_col, depth + 1, tri, CalculateColour);
        else
           RayTrianglesIntersection(p, refl, t, refl_col, depth + 1, tri, CalculateColour);
    }

    glm::vec3 ret = amb + diff + (refl_col * atten);
//    glm::vec3 ret = tri->v1.col + (refl_col * atten);

    return ret;
}
*/

vec3 random_unit_vector()
{
    while (true)
    {
        float randx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        randx *= 2.f;
        randx -= 1.f;

        float randy = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        randy *= 2.f;
        randy -= 1.f;

        float randz = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        randz *= 2.f;
        randz -= 1.f;

        vec3 v = vec3(randx, randy, randz);
        if (!(v.x == 0.f || v.y == 0.f || v.z == 0.f))
            return normalize(v);
    }
}


glm::vec3 CalculateColourFlat(triangle* tri, int depth, glm::vec3 p, glm::vec3 dir, payload *prd)
{
    prd->rayHitSky = 0;

    vec3 worldNrm = tri->v1.nor;

    vec3 offset = worldNrm; 
    offset *= (0.0001 * sign(dot(dir, worldNrm)));
    vec3 origin = p - offset;
    prd->rayOrigin = origin;

    // For a random diffuse bounce direction, we follow the approach of
    // Ray Tracing in One Weekend, and generate a random point on a sphere
    // of radius 1 centered at the normal. This uses the random_unit_vector
    // function from chapter 8.5:
    const float theta = 6.2831853 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);  // Random in [0, 2pi]
    const float u = 2.0 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 1.0;  // Random in [-1, 1]
    const float r = sqrt(1.0 - u * u);
    vec3 rayDir = worldNrm + vec3(r * cos(theta), r * sin(theta), u);
    // Then normalize the ray direction:
    rayDir = normalize(rayDir);

//    vec3 rayDir = worldNrm + random_unit_vector(); 

    prd->rayDir = rayDir;
    prd->hitNormal = worldNrm;

    vec3 colour = vec3(0.7f);

    //  if(worldNrm.y > .99999)
    //    prd.colour  = vec3(1, 0, 0);
    //  if(worldNrm.z > .9999999)
    //    prd.colour  = vec3(0, 1, 0);
    //  if(worldNrm.x < -.89999)
    //    prd.colour  = vec3(0, 0, 1);
    //  if(worldNrm.x > .9999999)
    //    prd.colour  = vec3(0, 1, 0);

//    if (tri->v1.nor.x > 0.99)
//        colour = vec3(1, 0, 0);
//    if (tri->v1.nor.x < -0.99)
//        colour = vec3(0, 1, 0);
//    if (tri->v1.nor.y > 0.99)
//        colour = vec3(0, 0, 1);
    if (tri->primitiveID == 0 || tri->primitiveID == 1)
        colour = vec3(1, 0, 0);
    if (tri->primitiveID == 8 || tri->primitiveID == 9)
        colour = vec3(0, 1, 0);
    if (tri->primitiveID == 6 || tri->primitiveID == 7)
        colour = vec3(0, 0, 1);

    return colour;
}


/*
glm::vec3 CalculateColourBVH(triangle* tri, int depth, glm::vec3 p, glm::vec3 dir)
{
    glm::vec3 refl_col(0.f);
    float t;

    if (tri->reflect && depth < max_recursion_depth)
    {
        glm::vec3 refl = glm::reflect(dir, tri->v1.nor);
        RayBVHIntersection(p, refl, &g_BVH, t, refl_col, depth + 1, tri);
    }

    glm::vec3 ret = tri->v1.col + (refl_col * atten);

    return ret;
}
*/
/*
void RayTrianglesIntersection(glm::vec3 o, glm::vec3 dir, float& t, glm::vec3 &io_col, int depth, triangle *from_tri, closest_hit p_hit)
{
    float closest_t = FLT_MAX;
    int closest_tc = -1;
    glm::vec3 closest_p;
    for (int tc = 0; tc < tris.size(); tc++)
    {
        triangle* tri = &tris[tc];

        if (from_tri != NULL && tri == from_tri)
            continue;

        glm::vec3 p = glm::vec3(0);
        float current_t = RayTriangleIntersection(o, dir, tri, p);

        if (current_t > 0 && current_t < f)
        {
            if (current_t < closest_t)
            {
                closest_t = current_t;
                closest_tc = tc;
                closest_p = p;
            }
        }
    }

    if (closest_tc >= 0)
    {
        t = closest_t;
        triangle* closest_tri = &tris[closest_tc];
        io_col = p_hit(closest_tri, depth, closest_p, dir, NULL);
        return;
    }

    io_col = bkgd;
    return;
}

void RayTraceTriangles()
{
    for (int pixel_y = 0; pixel_y < PIXEL_H; ++pixel_y)
    {
        float percf = (float)pixel_y / (float)PIXEL_H;
        int perci = percf * 100;
        std::clog << "\rScanlines done: " << perci << "%" << ' ' << std::flush;

        for (int pixel_x = 0; pixel_x < PIXEL_W; ++pixel_x)
        {

            glm::vec3 pp = GetPixelInViewSpace(pixel_x, pixel_y, PIXEL_W, PIXEL_H, l, r, t, b, n, f);
            glm::vec3 dir = normalize(pp);
            dir.y = -dir.y;
            dir.z = -dir.z;

            float t = 0;
            glm::vec3 col(0.f);            
            RayTrianglesIntersection(eye, dir, t, col, 0, NULL, CalculateColour);

//            if (t > 0 && t < f)
            {
                float pixel_r = col.x > 1 ? 255 : col.x * 255;
                float pixel_g = col.y > 1 ? 255 : col.y * 255;
                float pixel_b = col.z > 1 ? 255 : col.z * 255;

                pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 0] = pixel_r;
                pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 1] = pixel_g;
                pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 2] = pixel_b;

//                SDL_SetRenderDrawColor(renderer, pixel_r, pixel_g, pixel_b, 255);
//                SDL_RenderDrawPoint(renderer, pixel_x, pixel_y);
            }
        }
    }

    std::clog << "\rFinish rendering.           \n";

}
*/

void swap(float& l, float& r)
{
    float t = l;
    l = r;
    r = t;
}

bool RayAABBIntersection(glm::vec3 o, glm::vec3 dir, AABB* bb)
{
    float tmin = (bb->ftl.x - o.x) / dir.x;
    float tmax = (bb->bbr.x - o.x) / dir.x;

    if (tmin > tmax)
    {
        swap(tmin, tmax);
    }

    float tymin = (bb->ftl.y - o.y) / dir.y;
    float tymax = (bb->bbr.y - o.y) / dir.y;

    if (tymin > tymax)
    {
        swap(tymin, tymax);
    }

    if ((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    float tzmin = (bb->ftl.z - o.z) / dir.z;
    float tzmax = (bb->bbr.z - o.z) / dir.z;

    if (tzmin > tzmax) 
    {
        swap(tzmin, tzmax);
    }

    if ((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }

    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    return true;
}


void RayBVHIntersection(glm::vec3 o, glm::vec3 dir, 
    BVH_node* inBVH, 
    float &io_t, glm::vec3 &ioCol, 
    int depth, triangle* from_tri, 
    payload *io_prd)
{
    vec3 miss_col = skyColour(dir);

    glm::vec3 refl_col(0.f);

    bool i = RayAABBIntersection(o, dir, &inBVH->bb);

    if (i == FALSE)
    {
        ioCol = miss_col;
        io_prd->rayHitSky = 1;
        return;
    }

    //no child BVH so test the triangles
    if (inBVH->child_bvhs.size() <= 0)
    {
        glm::vec3 p = glm::vec3(0);
        glm::vec3 closest_p = glm::vec3(0);
        float t = FLT_MAX;
        float closest_t = FLT_MAX;
        triangle* closest_tri = NULL;
        payload closest_prd;

        for (auto tri : inBVH->tris)
        {
            if (from_tri != NULL && tri == from_tri)
                continue;
            
            t = RayTriangleIntersection(o, dir, tri, p);
            if (t < closest_t)
            {
                closest_t = t;
                closest_tri = tri;
                closest_p = p;
            }
        }

        if (closest_t == FLT_MAX)
        {
            ioCol = miss_col;
            io_prd->rayHitSky = 1;
            return;
        }

        io_t = closest_t;
        ioCol = CalculateColourFlat(closest_tri, depth, closest_p, dir, &closest_prd);
        io_prd->rayOrigin = closest_prd.rayOrigin;
        io_prd->rayDir = closest_prd.rayDir;
        io_prd->rayHitSky = closest_prd.rayHitSky;
        return;
    }



    float t = FLT_MAX;
    payload prd;
    float closest_t = FLT_MAX;
    glm::vec3 col(0.f);
    glm::vec3 closest_col(0.f);
    payload closest_prd;
    for (auto bvh : inBVH->child_bvhs)
    {
        RayBVHIntersection(o, dir, bvh, t, col, depth, from_tri, &prd);
        if (t < closest_t)
        {
            closest_t = t;
            closest_col = col;
            closest_prd.rayOrigin = prd.rayOrigin;
            closest_prd.rayDir = prd.rayDir;
            closest_prd.rayHitSky = prd.rayHitSky;
        }
    }

    if (closest_t == FLT_MAX)
    {
        ioCol = miss_col;
        io_prd->rayHitSky = 1;
        return;
    }

    io_t = closest_t;
    ioCol = closest_col;
    io_prd->rayOrigin = closest_prd.rayOrigin;
    io_prd->rayDir = closest_prd.rayDir;
    io_prd->rayHitSky = closest_prd.rayHitSky;
    return;


    ioCol = miss_col;
    return;
}

void PathTraceBVH()
{
    payload prd;
    for (int pixel_y = 0; pixel_y < PIXEL_H; ++pixel_y)
    {
        float percf = (float)pixel_y / (float)PIXEL_H;
        int perci = percf * 100;
        std::clog << "\rScanlines done: " << perci << "%" << ' ' << std::flush;

        for (int pixel_x = 0; pixel_x < PIXEL_W; ++pixel_x)
        {

            glm::vec3 summed_pixel_color(0, 0, 0);
            for (int sample = 0; sample < samples_per_pixel; sample++)
            {

                glm::vec3 pp = GetRandomPixelPosInViewSpace(pixel_x, pixel_y, PIXEL_W, PIXEL_H, l, r, t, b, n, f);

                vec3 origin = eye;

                glm::vec3 dir = normalize(pp);
                dir.y = -dir.y;
                dir.z = -dir.z;

                vec3 accumulatedRayColor = vec3(1.0);  // The amount of light that made it to the end of the current ray.

//                origin = prd.rayOrigin;
//                dir = prd.rayDir;

                for (int segments = 0; segments < tracedSegments; segments++)
                {

                    float t = -1;
                    glm::vec3 col(0.f);
                    RayBVHIntersection(origin, dir, &g_BVH, t, col, 0, NULL, &prd);

                    accumulatedRayColor *= col;

                    if (prd.rayHitSky == 1)
                    {
                        summed_pixel_color += accumulatedRayColor;
                        break;
                    }
                    else
                    {
                        origin = prd.rayOrigin;
                        dir = prd.rayDir;
                        prd.rayHitSky = 0;  // Will stop if light is hit
                    }
                }
            }

            vec3 average_pixel_colour = summed_pixel_color / float(samples_per_pixel);

            auto r = linear_to_gamma(average_pixel_colour.x);
            auto g = linear_to_gamma(average_pixel_colour.y);
            auto b = linear_to_gamma(average_pixel_colour.z);

            float pixel_r = r > 1 ? 255 : r * 255;
            float pixel_g = g > 1 ? 255 : g * 255;
            float pixel_b = b > 1 ? 255 : b * 255;

            pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 0] = pixel_r;
            pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 1] = pixel_g;
            pixelBuffer[(pixel_y * PIXEL_W * 3) + (pixel_x * 3) + 2] = pixel_b;


            //printf("%f\n", float(pixel_y * PIXEL_W + pixel_x) / float(PIXEL_H * PIXEL_W));
        }
    }
    std::clog << "\rFinish rendering.           \n";

}


void CounterEndAndPrint(LARGE_INTEGER StartingTime, LARGE_INTEGER *EndingTime, LARGE_INTEGER Frequency)
{
    QueryPerformanceCounter(EndingTime);
    
    LARGE_INTEGER ElapsedMicroseconds;
    ElapsedMicroseconds.QuadPart = (*EndingTime).QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout << (float)ElapsedMicroseconds.QuadPart / (float)1000000 << std::endl;
}

int main()
{
    srand((unsigned int)time(NULL));

    LARGE_INTEGER Frequency;
    QueryPerformanceFrequency(&Frequency);

//    const std::string MODEL_PATH = "objs/white_oak/white_oak.obj";
//    const std::string MODEL_PATH = "objs/bird/textured_quad.obj";
 //   const std::string MODEL_PATH = "objs/room/viking_room.obj";
        const std::string MODEL_PATH = "objs/cornell/Cornell.obj";
//    const std::string MODEL_PATH = "objs/sphere/sphere.obj";
//    const std::string MODEL_PATH = "objs/pokeball/pokeball.obj";

    obj_parse(MODEL_PATH.c_str(), &objs, .1f);

    tris = AssemblePrimitives(verts, n_verts);

//    AppendTriangles(&tris, objs);

    LARGE_INTEGER Construct_StartingTime, Construct_EndingTime;
    QueryPerformanceCounter(&Construct_StartingTime);

    if (use_bvh)
    {
        ConstructBVH(&g_BVH, tris, bvh_width);
    }

    CounterEndAndPrint(Construct_StartingTime, &Construct_EndingTime, Frequency);


    LARGE_INTEGER Render_StartingTime, Render_EndingTime;
    QueryPerformanceCounter(&Render_StartingTime);

    if (use_bvh)
    {
        PathTraceBVH();
    }
    else
    {
//        RayTraceTriangles();
    }

    CounterEndAndPrint(Render_StartingTime, &Render_EndingTime, Frequency);

    savebitmap("render.bmp", pixelBuffer, PIXEL_W, PIXEL_H);
   // savehdr("render.hdr", pixelBuffer, PIXEL_W, PIXEL_H);
    stbi_write_hdr("out.hdr", PIXEL_W, PIXEL_H, 3, pixelBuffer);

    return EXIT_SUCCESS;
}
