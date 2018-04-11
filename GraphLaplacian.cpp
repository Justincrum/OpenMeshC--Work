// GraphLaplacian.cpp : Defines the entry point for the console application.
//

#include <iostream>  //Including the required packages to run the code.
#include <vector>//Required package
#include <numeric>

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh> //Required package
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  //Required package
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh; //Typedef defines a simpler name for a data type.
int main(int argc, char **argv)  /* int main tells the program the output should be of integer type.
								 int argc tells the program to let me input c integer-valued arguments
								 from the command line.  char **argv tells stores the inputs as an
								 array.  If I do not want to run from command line, I can skip these?*/
{

	std::vector<std::vector<int>> Laplacian;  //Creates an initial storage location for the graph Laplacian.
	int degree = 0; //Starting with 0 vertices.
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
		MyMesh::Scalar              valence;  //Tells the iterations to also determine the number of neighbors.
		//unsigned int                i, N(atoi(argv[1])); //Best guess is this is taking in the input number of smoothing iterations we want.

		//for (i = 0; i < N; ++i) //Cycle through the number of smoothing iterations.
		//{
		//This following line is the iterator.
			for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) //Cycle through each vertex in the mesh.
			{
				++degree; //increase degree by 1 for each vertex in the mesh.
			}

			//Next loop reshapes the vector of vectors to be a square matrix of size degree x degree.
			for (int i = 0; i < degree; i++) {
				Laplacian[i].resize(degree);
			}

			for (int i = 0; i < degree; i++) {
				for (int j = 0; j < degree; j++) {
					Laplacian[i][j] = 0;
				}
			}

			for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) //Cycle through each vertex in the mesh.
			{
							   //The following line is the circulator that checks for neighbors of the original vertex.
				for (vv_it = mesh.vv_iter(*v_it); vv_it; ++vv_it) /* This loop counts the valence of each vertex. */
				{
					Laplacian[v_it][vv_it] = 1;
				}
			}
		
			std::vector<int> val; //Create a storage location for the valence of a vertex.
			//Compute the valence of each vertex.
			for (int i = 0; i < degree; i++) {
				int sum = std::accumulate(Laplacian[i].begin(), Laplacian[i].end(), 0);
				val[i] = sum;
				}
			//In the diagonal, store the negative of the valence of the vertex.
			for (int i = 0; i < degree; i++) {
				for (int j = 0; j < degree; j++) {
					if (i == j)
						Laplacian[i][j] = -val[i];
				}
			}

		//}
		return 0; //If everything runs correctly it returns 0 to symbolize no errors.
	}
