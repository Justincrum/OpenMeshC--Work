#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include "eigensolver.h"


using namespace std;

#define WIDTH 600
#define HEIGHT 600
#define PAIRWIDTH 0.02/20
#define EDGEWIDTH 0.015/20
#define EIGENVECTORNUMS 75


// debug, output time cost
//------------------------
#include <mach/mach_time.h>

double conversion_factor;

void Init() {
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    conversion_factor = (double)timebase.numer / (double)timebase.denom / 1e09;
}
//-----------------------

using namespace std;

std::vector<double> eigenvalues;
std::vector<Eigen::VectorXd> eigenvectors;
int current_eigen_index = 1;



class MyTraits : public OpenMesh::DefaultTraits
{
public:
    //typedef OpenMesh::Vec3d Point; // use double-values points
    //VertexAttributes(OpenMesh::Attributes::Status);
    //FaceAttributes(OpenMesh::Attributes::Status);
    //EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits>  MyPolyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;


MyMesh mesh;

OpenMesh::HPropHandleT<double> Angle;               // actually store cot, not the angle
OpenMesh::VPropHandleT<MyMesh::Normal> Lap_Bel;
OpenMesh::VPropHandleT<double> EigenFunc;




void computeAngle(MyMesh & mesh){
    for (auto h_it = mesh.halfedges_begin(); h_it!= mesh.halfedges_end(); ++h_it){
        MyMesh::HalfedgeHandle next = mesh.next_halfedge_handle(*h_it);
        MyMesh::FaceHandle fh = mesh.face_handle(*h_it);
        if (fh.is_valid()){
            MyMesh::Normal v1 = mesh.calc_edge_vector(*h_it).normalize();
            MyMesh::Normal v2 = mesh.calc_edge_vector(next).normalize();
            mesh.property(Angle,*h_it) = (-1.0 * v1).operator|(v2) / (v2 % v1).norm();
        }//else{ std::cout << "not valid" << std::endl;}
    }
}

void Laplace_Operation(MyMesh & mesh){
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
        for (auto vih_it = mesh.vih_iter(*v_it); vih_it.is_valid(); ++vih_it){
            MyMesh::HalfedgeHandle next = mesh.next_halfedge_handle(*vih_it);
            if (mesh.face_handle(*vih_it).is_valid()){
                MyMesh::Normal v_ij = mesh.calc_edge_vector(*vih_it);
                mesh.property(Lap_Bel,*v_it) += 0.5 * mesh.property(Angle,next) * v_ij;
            }
        }
        for (auto voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); ++voh_it){
            MyMesh::HalfedgeHandle next = mesh.next_halfedge_handle(*voh_it);
            if (mesh.face_handle(*voh_it).is_valid()){
                MyMesh::Normal v_ji = mesh.calc_edge_vector(*voh_it);
                mesh.property(Lap_Bel,*v_it) -= 0.5 * mesh.property(Angle,next) * v_ji;
            }
        }
    }
}

void Laplace_Matrix(Eigen::MatrixXd & matrix, MyMesh & mesh){
    for (auto h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it){

        MyMesh::FaceHandle fh = mesh.face_handle(*h_it);
        if (!fh.is_valid()) continue;

        MyMesh::VertexHandle to_vertex = mesh.to_vertex_handle(*h_it);
        MyMesh::VertexHandle from_vertex = mesh.from_vertex_handle(*h_it);
        int i = to_vertex.idx();
        int j = from_vertex.idx();
        matrix(i,j) -= 0.5;// * mesh.property(Angle,mesh.next_halfedge_handle(*h_it));
        matrix(j,i) -= 0.5;// * mesh.property(Angle,mesh.next_halfedge_handle(*h_it));
        matrix(i,i) += 0.5;// * mesh.property(Angle,mesh.next_halfedge_handle(*h_it));
        matrix(j,j) += 0.5;// * mesh.property(Angle,mesh.next_halfedge_handle(*h_it));
    }
}

void Laplace_Matrix(Eigen::SparseMatrix<double> & matrix, MyMesh & mesh){
    std::vector< Eigen::Triplet<double>  > tripletList;
    tripletList.reserve(mesh.n_edges()*2*4);

    for (auto h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it){

        MyMesh::FaceHandle fh = mesh.face_handle(*h_it);
        if (!fh.is_valid()) continue;

        MyMesh::VertexHandle to_vertex = mesh.to_vertex_handle(*h_it);
        MyMesh::VertexHandle from_vertex = mesh.from_vertex_handle(*h_it);
        int i = to_vertex.idx();
        int j = from_vertex.idx();
        tripletList.push_back( Eigen::Triplet<double> (i,j,-0.5 * mesh.property(Angle,mesh.next_halfedge_handle(*h_it))) );
        tripletList.push_back( Eigen::Triplet<double> (j,i,-0.5 * mesh.property(Angle,mesh.next_halfedge_handle(*h_it))) );
        tripletList.push_back( Eigen::Triplet<double> (i,i,0.5 * mesh.property(Angle,mesh.next_halfedge_handle(*h_it))) );
        tripletList.push_back( Eigen::Triplet<double> (j,j,0.5 * mesh.property(Angle,mesh.next_halfedge_handle(*h_it))) );
    }
    matrix.resize(mesh.n_vertices(),mesh.n_vertices());
    matrix.setFromTriplets(tripletList.begin(),tripletList.end());
}





int main(int argc, char ** argv){

    OpenMesh::IO::Options ropt, wopt;
    mesh.request_vertex_normals();
    //mesh.request_face_normals();
    ropt += OpenMesh::IO::Options::VertexNormal;

    if (!OpenMesh::IO::read_mesh(mesh, argv[1], ropt))
    {
        std::cerr << "read error\n";
        exit(1);
    }

    if (!ropt.check(OpenMesh::IO::Options::VertexNormal)){
        std::cout << "no normal information, compute normals" << std::endl;
        mesh.request_face_normals();
        mesh.update_normals();
    }else{
        mesh.request_face_normals();
        mesh.update_face_normals();
    }


    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(mesh.n_vertices(),mesh.n_vertices());
    mesh.add_property(Angle);
    mesh.add_property(Lap_Bel);
    mesh.add_property(EigenFunc);
    computeAngle(mesh);
    Laplace_Matrix(matrix,mesh);

    Eigen::SparseMatrix<double> sparse_matrix;
    Laplace_Matrix(sparse_matrix,mesh);

    //std::cout << sparse_matrix << std::endl;

    // debug , time cost
    // ------
    Init();
    uint64_t t1, t2;
    t1 = mach_absolute_time();
    // -----

    SparseEigenSolver solver(sparse_matrix,EIGENVECTORNUMS);
    solver.computeSpectrum();

    t2 = mach_absolute_time();
    double duration_s = (double)(t2 - t1) * conversion_factor;
    std::cout << "compute eigenfunctions, cost time " << duration_s << std::endl << flush;


    eigenvalues = solver.getEigenvalues();
    for (int i = 0; i < eigenvalues.size(); ++i)
    std::cout << "eigenvalue " << i << " "<<eigenvalues[i] << std::endl;
    eigenvectors = solver.getEigenvectors();
    //std::cout << eigenvectors[1] << std::endl;

    for (int i = 0; i < eigenvectors[current_eigen_index].size(); i++){
        MyMesh::VertexHandle vh(i);
        mesh.property(EigenFunc,vh) = eigenvectors[current_eigen_index](i);
    }


    return 0;
}

