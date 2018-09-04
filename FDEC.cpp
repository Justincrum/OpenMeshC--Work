//Justin Crum
//8/23/2018

/*FDEC is a code to implement attempt 1 of fractional discrete exterior calculus.
Note that on a theoretical level, we have no fractional generalized Stokes' theorem to use
in discretizing a fractional k-form. */



#include <utility>
#include <iostream>  
#include <vector>
#include <numeric>
#include <set>
#include<Eigen/Dense>




// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh> //Required package
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>



//Defining the meshes we will use.
//typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh; 
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

MyMesh mesh;

//Now we create some vertex iterators.
MyMesh::VertexIter   v_it, v_end(mesh.vertices_end());
MyMesh::VertexVertexIter vv_it;

//Function that will be used to generate random values between 1 and 20 for the
//0-form ZeroF.

unsigned int RandomValue(int x,int y)
{
	srand(y);

	return rand() % x + 1;
}

//This is the weight function W(J,alpha).  We can modify this as we find necessary.
double Weight(double x, double y)
{
	return (1 - y) / x;
}

//In the following three functions, nv will be the number of vertices.
//This is what we will use to choose out vertices that are a minimum distance.
int minDistance(std::vector<int> dist, std::vector<bool> sptSet)
{
	int min = INT_MAX, min_index;

	for (int v = 0; v < mesh.n_vertices(); v++)
	{
		if (sptSet[v] == false && dist[v] <= min)
		{
			min = dist[v], min_index = v;
		}
	}

	return min_index;
}


//Use this function to check that it is working correctly.
void printSolution(std::vector<int> dist, int n)
{
	printf("Vertex Distance from Source \n");
	for (int i = 0; i < mesh.n_vertices(); i++) 
	{
		printf("%d    %d\n", i, dist[i]);
	}
}
//Not sure how to make use of this.
void Graph()
{
	std::vector<std::vector<int>> graph;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) //Cycle through each vertex in the mesh.
	{
		//The following line is the circulator that checks for neighbors of the original vertex.
		for (vv_it = mesh.vv_iter(*v_it); vv_it; ++vv_it) //This loop inputs a value of 1 at each neighbor.
		{
			int m = v_it->idx();
			int n = vv_it->idx();
			graph[m][n] = 1;
		}
	}


}

void Dijkstra(std::vector<std::vector<int>> graph, int src)
{
	std::vector<int> dist;
	std::vector<bool> sptSet;

	
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		dist[i] = INT_MAX, sptSet[i] = false;
	}

	dist[src] = 0;

	for (int count = 0; count < mesh.n_vertices()-1; count++)
	{
		int u = minDistance(dist, sptSet);

		sptSet[u] = true;

		for (int v = 0; v < mesh.n_vertices(); v++)
		{
			if (!sptSet[v] && graph[u][v] == 1 && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v])
			{
				dist[v] = dist[u] + graph[u][v];
			}
		}

		printSolution(dist, mesh.n_vertices());
	}
}


int main(int argc, char **argv) 

{

	//Checking for correct inputs.
	if (argc != 3) 
	{			  
		           //Want three inputs: Input 0 = code name
				   //					Input 1 = Mesh file
				   //					Input 2 = Alpha value

		std::cerr << "Usage: " << argv[0] << "#infile\n";

		return 1;

	}

	//Now we read in the mesh.

	if (!OpenMesh::IO::read_mesh(mesh, argv[1]))
	{
		std::cerr << "Error:  Cannot read mesh from: " << argv[1] << std::endl;
		return 1;

	}

	//Create a custom vertex property that stores a value for each vertex, overall
	//this is our 0-form, which we call ZeroF.
	OpenMesh::VPropHandleT<MyMesh::Scalar> ZeroF;
	mesh.add_property(ZeroF);

	unsigned int numberVertices = mesh.n_vertices();

	//Now we create some vertex iterators.
	//MyMesh::VertexIter   v_it, v_end(mesh.vertices_end());
	//MyMesh::VertexVertexIter vv_it;

	//Stores a random value between 1 and 20 at each vertex.

	////////Why doesn't this part work?  Gives error about a circulator in circulatorst.hh line 457
	//for (v_it = mesh.vertices_begin(); v_it != v_end; ++vv_it)
	//{
	//	mesh.property(ZeroF, v_it) = RandomValue(20);
	//}


	//Create the adjacency matrix.  Need this for Dijkstra's algorithm.
	std::vector<std::vector<int>> graph;

	MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());  //This is what allows me to iterate through each vertex (?)
	MyMesh::VertexVertexIter    vv_it;  /*VertexVertexIter I think will sub-iterate through at each step of VertexIter.
										So this line lets me access the neighbors of each vertex in the in mesh.*/

	int nv;
	//Count the number of vertices on the mesh.
	nv = mesh.n_vertices();

	//The next 4 loops are used to create the graph adjacency matrix.  
	Eigen::MatrixXf Laplacian(nv, nv);
	for (int i = 0; i < nv; i++) {
		for (int j = 0; j < nv; j++) {
			Laplacian(i, j) = 0;
		}
	}

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) //Cycle through each vertex in the mesh.
	{
		//The following line is the circulator that checks for neighbors of the original vertex.
		for (vv_it = mesh.vv_iter(*v_it); vv_it; ++vv_it) //This loop inputs a value of 1 at each neighbor.
		{
			Laplacian(v_it->idx(), vv_it->idx()) = 1;
		}
	}

	graph.resize(nv);
	for (int i = 0; i < nv; i++)
	{
		graph[i].resize(nv);
	}

	for (int i = 0; i < nv; i++)
	{
		for (int j = 0; j < nv; j++)
		{
			graph[i][j] = Laplacian(i, j);
		}
	}

	std::vector<std::vector<int>> DistanceArray;

	DistanceArray.resize(nv);
	for (int j = 0; j < nv; j++)
	{
		DistanceArray[j].resize(nv);
	}

	//The following loop creates a square matrix that we can use to read off distances between
	//two vertices.  
	
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{

		//These will store the points v_start and v_end.
		//MyMesh::VertexHandle endV = mesh.to_vertex_handle(mesh.halfedge_handle(e_it, 0));
		//MyMesh::VertexHandle startV = mesh.from_vertex_handle(mesh.halfedge_handle(e_it, 0));

		//int locEndV = endV.idx();
		//int locStartV = startV.idx();

		int locV = v_it ->idx();

		std::vector<int> dist;
		dist.resize(nv);

		std::vector<bool> sptSet;
		sptSet.resize(nv);

		int src = locV;

		for (int i = 0; i < nv; i++)
		{
			dist[i] = INT_MAX, sptSet[i] = false;
		}

		dist[src] = 0;

		for (int count = 0; count < mesh.n_vertices() - 1; count++)
		{
			int u = minDistance(dist, sptSet);

			sptSet[u] = true;

			for (int v = 0; v < mesh.n_vertices(); v++)
			{
				if (!sptSet[v] && graph[u][v] == 1 && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v])
				{
					dist[v] = dist[u] + graph[u][v];
				}
			}

			//This loop sets up the matrix of distances.  Recall that the way the matrix is used is as
			//follows:
			//to get d(v_j, v_0) or d(v_0, v_j), use zeroth row.
			//to get anything else, you can choose d(v_i, v_end) = a_(i,idx(v_end)) as the entry.
			//The zeroth column keeps track of which vertex we are thinking about.
			for (int i = 0; i < nv; i++)
			{
				DistanceArray[locV][i] = dist[i];
				DistanceArray[i][0] = i;

			}
		}

	}
	//Now that we have the array of distances, we can implement our own schemes!



	//Now we want to implement idea for fractional discrete exterior derivative
	//from vertex values to edge values.  

	//The list of edge values defines the discrete (1,alpha)-form d^alpha ZeroF,
	//where ZeroF is our discrete 0-form.

	////Want to make a 1-d array of ZeroF's value at each point.
	std::vector<unsigned int> ZeroFList;

	//seed number.
	int j = 0;

	//Assign values to each vertex.
	for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{

		ZeroFList.push_back(RandomValue(20, j));
		++j;
	}

	//Print ZeroF values at each vertex for checking purposes.
	for (int i = 0; i < numberVertices; ++i)
	{
		std::cout << ZeroFList[i] << std::endl;
	}

	//
	double alpha;
	alpha = strtod(argv[2], NULL);

	std::vector<double> DZeroF;

	//Now that we have a 0-form ZeroF and the graph distance of any two pairs of vertices, we can
	//implement our own algorithm to try to understand what it is doing.

	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		//Create indices for the beginning and end of the edge.
		MyMesh::VertexHandle endV = mesh.to_vertex_handle(mesh.halfedge_handle(e_it, 0));
		MyMesh::VertexHandle startV = mesh.from_vertex_handle(mesh.halfedge_handle(e_it, 0));
		int locEndV = endV.idx();
		int locStartV = startV.idx();

		//Now we want to create two sets, say E and S.  The set E is a list of all vertices
		//closer to endV, and the set S is the set of all vertices closer to startV.
		//Any vertex that is equidistant does not contribute.  Both sets E and S will only contain
		//the index of a vertex.  

		std::vector<int> S;
		std::vector<int> E;

		double sumE = 0;
		double sumS = 0;

		for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			int j = v_it->idx();

			//Notice that we do not partition the whole mesh--equidistant vertices do not contribute.
			if (DistanceArray[locStartV][j] < DistanceArray[locEndV][j] && DistanceArray[locStartV][j] > 0)
			{
				S.push_back(v_it->idx());
			}	

			if (DistanceArray[locStartV][j] > DistanceArray[locEndV][j] && DistanceArray[locEndV][j] > 0)
			{
				E.push_back(v_it->idx());
			}

		}
		//Now that we have the sets E and S, we compute:

		while (!E.empty())
		{
			if (E[0] == 0)
			{
				sumE = sumE + (ZeroFList[E[0]] * Weight(DistanceArray[E[0]][locEndV], alpha));
				E.erase(E.begin());
			} 
			else
			{
				sumE = sumE + (ZeroFList[E[0]] * Weight(DistanceArray[locEndV][E[0]], alpha));
				E.erase(E.begin());
			}

		}

		while (!S.empty())
		{
			if (S[0] == 0)
			{
				sumS = sumS + (ZeroFList[0] * Weight(DistanceArray[S[0]][locStartV], alpha));
				S.erase(S.begin());
			} 
			else
			{
				sumS = sumS + (ZeroFList[S[0]] * Weight(DistanceArray[locStartV][S[0]], alpha));
				S.erase(S.begin());
			}

			
		}

		sumE = sumE + ZeroFList[locEndV];
		sumS = sumS + ZeroFList[locStartV];

		//This populates the values at the edges based off the algorithm for the
		//fractional discrete alpha-derivative.
		DZeroF.push_back(sumE - sumS);
		std::cout << sumE-sumS << std::endl;




	}




} 