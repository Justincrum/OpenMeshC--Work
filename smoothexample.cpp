#include <iostream>  //Including the required packages to run the code.
#include <vector>//Required package
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh> //Required package
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  //Required package
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh; //Typedef defines a simpler name for a data type.
int main(int argc, char **argv)  /* int main tells the program the output should be of integer type.
								 int argc tells the program to let me input c integer-valued arguments
								 from the command line.  char **argv tells stores the inputs as an 
								 array.  If I do not want to run from command line, I can skip these?*/
{
    MyMesh  mesh; //I think this creates a base data type that we call mesh.  This is what we will add vertices and edges to.
    // check command line options
    if (argc != 4) //Checking if the number of arguments being input is not equal to 4.
    {
        std::cerr << "Usage:  " << argv[0] << " #iterations infile outfile\n"; //std is a built in namespace of functions in c++
		//std::cerr calls function cerr from std.  cerr stands for console error and is used to send an error message.
        return 1;  //1 represents an error occured.
    }
    // read mesh from stdin
    if ( ! OpenMesh::IO::read_mesh(mesh, argv[2]) )  //Making sure the program has a mesh to read (?)
    {
        std::cerr << "Error: Cannot read mesh from " << argv[2] << std::endl;  //If the above if statement fails to find a mesh, return an error.
        return 1; //1 symbolizes error.
    }
    // this vertex property stores the computed centers of gravity
    OpenMesh::VPropHandleT<MyMesh::Point> cogs;  /*This assigns a handle to the points of the mesh--cogs.*/
    mesh.add_property(cogs); /*property(...) allows me to add a custom property to the mesh.  mesh.add_property(cogs) tells open mesh that 
							 I want to be able to assign a cog value to each vertex.*/
    // smoothing mesh argv[1] times
    MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());  //This is what allows me to iterate through each vertex (?)
    MyMesh::VertexVertexIter    vv_it;  /*VertexVertexIter I think will sub-iterate through at each step of VertexIter. 
										So this line lets me access the neighbors of each vertex in the in mesh.*/
    MyMesh::Point               cog;  /*Tells the iterations to compute the center of gravity of the initial vertex based on
									  its neighboring vertices. */
    MyMesh::Scalar              valence;  //Tells the iterations to also determine the number of neighbors.
    unsigned int                i, N(atoi(argv[1])); //Best guess is this is taking in the input number of smoothing iterations we want.

    for (i=0; i < N; ++i) //Cycle through the number of smoothing iterations.
    {
        for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it) //Cycle through each vertex in the mesh.
        {
            mesh.property(cogs,*v_it).vectorize(0.0f); //Setting a storage spot for the center of gravity?
            valence = 0.0; //Starting valence at 0.

            for (vv_it=mesh.vv_iter( *v_it ); vv_it; ++vv_it) /* This loop does two things:  1. It adds together the position of each of
															  the neighbors of a vertex.  2. It increases the valence by 1 for each
															  neighbor the vertex has.*/
            {
                mesh.property(cogs,*v_it) += mesh.point( *vv_it );  /*+=mesh.point(..) takes the sum of the previous neighbors and adds
																	the next neighbor, each time as coordinates.*/
                ++valence; //++valence increases valence by 1 for each time we pass through the loop.
            }
            mesh.property(cogs,*v_it) /= valence;  /*Computes the actual center of gravity by taking the sum of the
												   coordinates of the neighbors and dividing by the computed valence in the inner loop. */
        }

        for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it) //Looping through vertices.
            if ( !mesh.is_boundary( *v_it ) )  //Checking if a mesh point is on the boundary.  
                mesh.set_point( *v_it, mesh.property(cogs,*v_it) ); /*Putting each point that is not on the boundary to the center of gravity?
																	This way if a point is on the boundary, the shape doesn't lose volume. */
    }
    // write mesh to stdout
    if ( ! OpenMesh::IO::write_mesh(mesh, argv[3]) )
    {
        std::cerr << "Error: cannot write mesh to " << argv[3] << std::endl;  //Checking if there is an appropriate location to store the output.
        return 1;  
    }
    return 0; //If everything runs correctly it returns 0 to symbolize no errors.

    MyMesh::FaceIter faceIter;  //This is iterating over face of mesh.
    faceIter = mesh.faces_begin(); //Start at face 1?
    MyMesh::FaceVertexIter fv = mesh.fv_iter(*faceIter);  //Reconfigure the faces based off how the mesh was changed by moving vertices?

}