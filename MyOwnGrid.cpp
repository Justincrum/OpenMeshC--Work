// MyOwnGrid.cpp : Defines the entry point for the console application.
//

#include <iostream>

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

using namespace std;
// ----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
// ----------------------------------------------------------------------------
// Build a n by m grid.

int main()
{
	MyMesh mesh;

	// generate vertices

	int* a = NULL; //Creates a pointer that will be filled in later.

	int n, m; //Creates spots for an input number of width vertices and length vertices.


	cin >> n; //Take in input values for number of y-values in grid. (rows)
	cin >> m; //Take in input values for number of x-values in grid. (columns)

	vector <MyMesh::VertexHandle> vhandle;


	for (int x = 0; x < m; x++) {
		for (int y = 0; y < n; y++) {
			vhandle.push_back(mesh.add_vertex(MyMesh::Point(x, y, 0)));
		}
		
	}
	

	// generate (quadrilateral) faces
	vector<MyMesh::VertexHandle>  face_vhandles;
	int iter = 0;
	
	for (int q = 0; q < m-1 ; q++) {
		for (int z = 0; z < n-1; z++) {
			face_vhandles.clear();
			face_vhandles.push_back(vhandle[iter + q]);
			face_vhandles.push_back(vhandle[iter + q + 1]);
			face_vhandles.push_back(vhandle[iter + q + 1 + n]);
			face_vhandles.push_back(vhandle[iter + q + n]);
			mesh.add_face(face_vhandles);
			iter++;
		}
		
	}
	


	// write mesh to gridnm.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "gridnm.obj"))
		{
			std::cerr << "Cannot write mesh to file 'gridnm.obj'" << std::endl;
			return 1;
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
		return 1;
	}
	return 0;
}




