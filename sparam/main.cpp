/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: chenh
 *
 * Created on July 30, 2018, 1:17 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <metis.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <fstream>

using namespace std;

/*
 * 
 */

typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
typedef OpenMesh::Vec3f Vec3D;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: \n"
                << argv[0] << "<input_mesh> <graph partition> \n\n";

        exit(1);
    }

    const bool OUTPUT_MESH_GRAPH = false;

    string filepath = std::string(argv[1]);

    // Reading mesh
    Mesh mesh_;
    mesh_.request_vertex_status();
    mesh_.request_vertex_normals();
    mesh_.request_face_normals();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    if (!OpenMesh::IO::read_mesh(mesh_, filepath)) {
        std::cerr << "OpenMesh failed to read mesh: "
                << argv[1] << std::endl;
        exit(1);
    }
    unsigned int n_vert = mesh_.n_vertices();
    unsigned int n_edge = mesh_.n_edges();


    // (optional) convert a mesh to graph format for metis partitioning
    if (OUTPUT_MESH_GRAPH) {

        Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
        Mesh::VertexVertexIter vv_it;

        ofstream graphstream;
        graphstream.open("temp_graph.txt");
        graphstream << n_vert << " " << n_edge << std::endl;

        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            for (vv_it = mesh_.vv_iter(v_it.handle()); vv_it.is_valid(); ++vv_it) {
                auto idx_n = vv_it->idx() + 1;
                graphstream << idx_n << " ";
            }
            graphstream << std::endl;
        }

        graphstream.close();
    }

    // read a partitioned file to find the list of vertex separators
    std::set<int> vert_separators_;
    std::set<int> vert_one_;
    std::set<int> vert_two_;    
    std::ifstream ifs(std::string(argv[2]), std::ifstream::in);
    
    if(ifs.good()){
        string line;
        int cidx = 0;
        while(std::getline(ifs, line))
        {
            int pidx = std::stoi(line);
            if(pidx == 0)
                vert_one_.insert(cidx);
            else
                vert_two_.insert(cidx);
            cidx++;
        }
        ifs.close(); 
    }
    else {
        std::cerr << "Failed to read partition files: "
                << argv[2] << std::endl;
        exit(1);
    }  

    if(vert_one_.size() + vert_two_.size() != n_vert)
    {
        std::cerr << "Bad partition: combined number of vertices should be " << n_vert
                << ", but it was read"<< vert_one_.size() + vert_two_.size() << std::endl;
        exit(1);
    }
    
    // find vertex separators
    for(auto idxit = vert_one_.begin(); idxit != vert_one_.end(); )
    {
        bool is_separator = false;
        OpenMesh::VertexHandle vh(*idxit);
        Mesh::VertexVertexIter vv_it;
        for (vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                auto idx_n = vv_it->idx();
                if(vert_two_.find(idx_n) != vert_two_.end())
                {
                    vert_separators_.insert(*idxit);
                    idxit = vert_one_.erase(idxit);
                    is_separator = true;
                    break;
                }               
        }
        
        if(!is_separator)
            ++idxit;
    }

    return 0;
}

