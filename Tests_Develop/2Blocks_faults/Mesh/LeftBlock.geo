Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box 
Box(1) = {0,0,0,2.5,10,15};
Physical Volume("Left", 1) = {1};
Physical Surface("Left_bound",1) = {1};
Physical Surface("Left_interf",2) = {2};
Physical Surface("Left_bottom",3) = {5};
Transfinite Surface {1,2,3,4,5,6};


//Recombine Surface {:};
MeshSize {:} = 0.75;
Mesh 3;
Save "LeftBlock.msh";
