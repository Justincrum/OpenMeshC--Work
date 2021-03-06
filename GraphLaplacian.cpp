// GraphLaplacian.cpp : Defines the entry point for the console application.
//

#include <utility>
#include <iostream>  //Including the required packages to run the code.
#include <vector>//Required package
#include <numeric>
#include <set>
#include<Eigen/Dense>


// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh> //Required package
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  //Required package

typedef std::pair <int, int> IntPair;
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh; //Typedef defines a simpler name for a data type.


int main(int argc, char **argv)  /* int main tells the program the output should be of integer type.
								 int argc tells the program to let me input c integer-valued arguments
								 from the command line.  char **argv tells stores the inputs as an
								 array.  If I do not want to run from command line, I can skip these?*/
{

	std::vector<std::vector<int>> Lap;  //Creates an initial storage location for the graph Laplacian.

		MyMesh  mesh; //I think this creates a base data type that we call mesh.  This is what we will add vertices and edges to.
					  // check command line options
		if (argc != 2) //Checking if the number of arguments being input is not equal to 4.
		{
			std::cerr << "Usage:  " << argv[0] << " #infile\n"; //std is a built in namespace of functions in c++
																				   //std::cerr calls function cerr from std.  cerr stands for console error and is used to send an error message.
			return 1;  //1 represents an error occured.
		}
		// read mesh from stdin
		if (!OpenMesh::IO::read_mesh(mesh, argv[1]))  //Making sure the program has a mesh to read (?)
		{
			std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;  //If the above if statement fails to find a mesh, return an error.
			return 1; //1 symbolizes error.
		}
								 //Checking neighbors of each vertex.
		MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());  //This is what allows me to iterate through each vertex (?)
		MyMesh::VertexVertexIter    vv_it;  /*VertexVertexIter I think will sub-iterate through at each step of VertexIter.
											So this line lets me access the neighbors of each vertex in the in mesh.*/

		int nv;
		//Count the number of vertices on the mesh.
		nv = mesh.n_vertices();

		//Initialize the Laplacian matrix with all 0 entries.

		Eigen::MatrixXd Laplacian(nv, nv);
		for (int i = 0; i < nv; i++) {
			for (int j = 0; j < nv; j++) {
				Laplacian(i, j) = 0;
			}
		}

		//Put 1's in each spot that indicates an edge between two vertices.

			for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) //Cycle through each vertex in the mesh.
			{
							   //The following line is the circulator that checks for neighbors of the original vertex.
				for (vv_it = mesh.vv_iter(*v_it); vv_it ; ++vv_it) //This loop inputs a value of 1 at each neighbor.
				{
					Laplacian(v_it->idx(), vv_it->idx()) = 1;
				}
			}

			std::vector<int> val; //Create a storage location for the valence of a vertex.
			val.resize(nv); //Put nv spots into the vector of valences.

			//Compute the valence of each vertex.
			for (int i = 0; i < nv; i++) {
				int sum = 0;
				for (int j = 0; j < nv; j++) {
					sum += Laplacian(i, j);
				}
				val[i] = sum;
				
				}
			//In the diagonal, store the negative of the valence of the vertex.
			for (int i = 0; i < nv; i++) {
				for (int j = 0; j < nv; j++) {
					if (i == j)
						Laplacian(i, j) = -val[i];

				}
			}


			std::cout << Laplacian << std::endl;
			
			std::cin.get();

		return 0; //If everything runs correctly it returns 0 to symbolize no errors.
	}
