#include "vcl/vcl.hpp"
#include "forest.hpp"
#include <iostream>

using namespace vcl;

static buffer<uint3> connectivity_grid(size_t Nu, size_t Nv)
{
    buffer<uint3> connectivity;
    for(size_t ku=0; ku<Nu-1; ++ku) {
        for(size_t kv=0; kv<Nv-1; ++kv) {
            unsigned int k00 = static_cast<unsigned int>(kv   + Nv* ku);
            unsigned int k10 = static_cast<unsigned int>(kv+1 + Nv* ku);
            unsigned int k01 = static_cast<unsigned int>(kv   + Nv*(ku+1));
            unsigned int k11 = static_cast<unsigned int>(kv+1 + Nv*(ku+1));

            connectivity.push_back(uint3{k00, k10, k11});
            connectivity.push_back(uint3{k00, k11, k01});
        }
    }
    return connectivity;
}

mesh create_terrain()

    {
    // Number of samples of the terrain is N x N
    const size_t N = 30;
    vec3 p00 = {-50,-0.5,-50};
    vec3 p10 = {50,-0.5,-50};
    vec3 p11 = {50,-0.5,50};
    vec3 p01 = {-50,-0.5,50};

    mesh terrain; // temporary terrain storage (CPU only)
    //terrain.position.resize(N*N);
    //terrain.color.resize(N*N);
    //terrain.texture.resize(N*N);

    // Fill terrain geometry
    for(size_t ku=0; ku<N; ++ku)
    {
        for(size_t kv=0; kv<N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku/(N-1.0f);
            const float v = kv/(N-1.0f);
            vec3 p = (1-u)*(1-v)*p00 + u*(1-v)*p10 + u*v*p11 + (1-u)*v*p01;
            vec3 const dpdu = (1-u)*(-p00+p01)+u*(-p10+p11);
            vec3 const dpdv = (1-v)*(-p00+p10)+v*( p11-p01);
            vec3 const n = normalize(cross( dpdv, dpdu));
            vec2 const uv = {u,v};

            const float scaling = 1.5f;
            const int octave = 2.0f;
            const float persistency = 3.0f;
            const float height = -0.5f;

            // Evaluate Perlin noise
            const float noise = noise_perlin(vec2(u,v), octave, persistency, 2.0f);

        
            // 3D vertex coordinates
            //float y = std::exp(-(u*u+v*v));
            float y = noise*height;
            p[1] = y;

        
            // Compute coordinates
            terrain.position.push_back(p);
            terrain.normal.push_back(n);
            terrain.uv.push_back(uv);
            
        }
    }



    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    terrain.connectivity = connectivity_grid(N,N);
    terrain.fill_empty_field();
    terrain.flip_connectivity();

    return terrain;
}

void generate_tree_positions(buffer<vec2>& positions){
    for (int i = 0; i<10; i++){
        float x = rand_interval(-50.0f,-25.0f);
        float z = rand_interval(-50.0f, 50.0f);
        positions.push_back({x,z});
    }

    for (int i = 0; i<10; i++){
        float x = rand_interval(25.0f,50.0f);
        float z = rand_interval(-50.0f, 50.0f);
        positions.push_back({x,z});
    }


}
