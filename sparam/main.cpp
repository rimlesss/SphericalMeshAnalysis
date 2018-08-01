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
#include <math.h>
#include <metis.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <Eigen/Sparse>
#include <QApplication>
#include <QtGui>
#include <QWidget>

using namespace std;

/*
 * 
 */

typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
typedef OpenMesh::Vec3f Vec3D;
typedef Eigen::Triplet<float> Trip;
typedef std::vector<Trip> Triplets;

Mesh mesh_;
std::set<int> vert_separators_;
std::set<int> vert_one_;
std::set<int> vert_two_;
std::map<int, int> index_map_;
std::map<int, float> bducoords_;
std::map<int, float> bdvcoords_;

class DrawWidget : public QWidget {
public:
    DrawWidget(const std::vector<float> pts, const std::vector<float> lines)
    {
        _pts = pts;
        _lines = lines;
    }
protected:
    void paintEvent(QPaintEvent *event) {
        //create a QPainter and pass a pointer to the device.
        //A paint device can be a QWidget, a QPixmap or a QImage
        QPainter painter(this);
        
        //create a black pen that has solid line
        //and the width is 2.
        QPen myPen(Qt::black, 2, Qt::SolidLine);
        painter.setPen(myPen);
        //draw the lines        
        painter.drawLine(100, 100, 100, 1);
        
        //draw the points
        myPen.setColor(Qt::red);
        painter.drawPoint(110, 110);
    }
private:
    std::vector<float> _pts;
    std::vector<float> _lines;
};

/*
 * Map the vertex separators uniformly to the boundary of a unit disk  
 * input: index i of the vertex in the separator list 
 * THIS METHOD SEEMS PROBLEMATIC!!?
 */
std::pair<float, float> getbduv(int idx) {
    auto n = (float) vert_separators_.size();
    float theta = (float) idx / (float) n * 2.0f * M_PI;
    return pair<float, float>(cos(theta), sin(theta));
}

void compute_valence(std::map<int, int> &valm, const set<int> include, const set<int> exclude) {
    for (auto idxit = include.begin(); idxit != include.end(); ++idxit) {
        Mesh::VertexVertexIter vv_it;
        OpenMesh::VertexHandle vh(*idxit);
        int valence = 0;
        for (vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
            if (exclude.find(vv_it->idx()) == exclude.end()) {
                valence++;
            }
        }
        valm.insert(pair<int, int>(*idxit, valence));
    }
}

/** 
 * Map a partition to a unit disk with barycentric method. The list of vertex separators are used as common boundaries.
 * Naive selection for edge weight: Dij = 1 (Tutte 1960)
 **/
void map_to_disk(const set<int> include, const set<int> exclude) {
    // compose valence map
    std::map<int, int> val_map;
    compute_valence(val_map, include, exclude);
    compute_valence(val_map, vert_separators_, exclude);

    auto plength = include.size();

    // prepare for Linear System
    auto eigen_A = Eigen::SparseMatrix<float>(plength, plength);
    eigen_A.setZero();

    // compose matrix A
    Triplets tflist;
    for (int i = 0; i < plength; i++)
        tflist.push_back(Trip(i, i, 1));

    for (auto idxit = include.begin(); idxit != include.end(); ++idxit) {
        Mesh::VertexVertexIter vv_it;
        OpenMesh::VertexHandle vh(*idxit);
        int idxi = index_map_.at(*idxit);
        float invalence = 1.0f / (float) val_map.at(*idxit);
        for (vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
            if (include.find(vv_it->idx()) != include.end()) {
                int idxj = index_map_.at(vv_it->idx());
                tflist.push_back(Trip(idxi, idxj, invalence));
            }
        }
    }
    eigen_A.setFromTriplets(tflist.begin(), tflist.end());
    delete &tflist;

    // compose vector U,V from boundary points
    Eigen::VectorXf eigen_u(plength);
    Eigen::VectorXf eigen_v(plength);

    for (auto idxit = include.begin(); idxit != include.end(); ++idxit) {
        Mesh::VertexVertexIter vv_it;
        OpenMesh::VertexHandle vh(*idxit);
        int idxi = index_map_.at(*idxit);
        float invalence = 1.0f / (float) val_map.at(*idxit);
        float u = 0.0f;
        float v = 0.0f;
        for (vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
            if (vert_separators_.find(vv_it->idx()) != vert_separators_.end()) {
                int idxj = index_map_.at(vv_it->idx());
                auto sepuv = getbduv(idxj);
                u += invalence * sepuv.first;
                v += invalence * sepuv.second;
            }
        }
        eigen_u[idxi] = u;
        eigen_v[idxi] = v;
    }

}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: \n"
                << argv[0] << "<input_mesh> <graph partition> \n\n";

        exit(1);
    }
    const bool OUTPUT_MESH_GRAPH = false;

    string filepath = std::string(argv[1]);

    // Reading mesh
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

    std::ifstream ifs(std::string(argv[2]), std::ifstream::in);


    if (ifs.good()) {
        string line;
        int cidx = 0;
        int vtidx = 0;
        while (std::getline(ifs, line)) {
            int pidx = std::stoi(line);
            if (pidx == 0)
                vert_one_.insert(cidx);
            else {
                vert_two_.insert(cidx);
                index_map_.insert(pair<int, int>(cidx, vtidx));
                vtidx++;
            }
            cidx++;
        }
        ifs.close();
    } else {
        std::cerr << "Failed to read partition files: "
                << argv[2] << std::endl;
        exit(1);
    }


    if (vert_one_.size() + vert_two_.size() != n_vert) {
        std::cerr << "Bad partition: combined number of vertices should be " << n_vert
                << ", but it was read" << vert_one_.size() + vert_two_.size() << std::endl;
        exit(1);
    }

    // find vertex separators
    int sep_counter = 0;
    int vfidx_counter = 0;
    for (auto idxit = vert_one_.begin(); idxit != vert_one_.end();) {
        bool is_separator = false;
        OpenMesh::VertexHandle vh(*idxit);
        Mesh::VertexVertexIter vv_it;
        for (vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
            auto idx_n = vv_it->idx();
            if (vert_two_.find(idx_n) != vert_two_.end()) {
                vert_separators_.insert(*idxit);
                index_map_.insert(pair<int, int>(*idxit, sep_counter));
                idxit = vert_one_.erase(idxit);
                is_separator = true;
                sep_counter++;
                break;
            }
        }

        if (!is_separator) {
            index_map_.insert(pair<int, int>(*idxit, vfidx_counter));
            idxit++;
            vfidx_counter++;
        }
    }

    // TODO: Verify that the separator vertices can be mapped well to a unit disk boundary




    /** 
     * Map the two partitions to a unit disk with barycentric method. 
     * The list of vertex separators are used as common boundaries 
     **/

    // run Qt Application
    QApplication app(argc,argv);
    DrawWidget qwidget;
    qwidget.show();
    return app.exec();
}

