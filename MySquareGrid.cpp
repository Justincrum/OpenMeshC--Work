#include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// ----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
// ----------------------------------------------------------------------------
// Build a simple cube and write it to std::cout

int grid_main()
{
	MyMesh mesh;
	// generate vertices
	MyMesh::VertexHandle vhandle[9];
	vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, -1, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(-1, 0, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(0, -1, 0));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(0, 1, 0));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(1, -1, 0));
	vhandle[7] = mesh.add_vertex(MyMesh::Point(1, 0, 0));
	vhandle[8] = mesh.add_vertex(MyMesh::Point(1, 1, 0));

	// generate (quadrilateral) faces
	std::vector<MyMesh::VertexHandle>  face_vhandles;

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[8]);
	face_vhandles.push_back(vhandle[7]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[7]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);

	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "output.off"))
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
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
	

