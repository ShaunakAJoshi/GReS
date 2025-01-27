close all;
clear;

%diary('ElasticVSdp.txt')
% simple mechanical model 
% A cube of size 1x1x1m is fixed in the bottom face and a load is applied
% at the top face. Both a load in vertical and horizontal direction is
% considered. Note that Drucker-Prager plasticity will not activate if only
% a vertical load is considered.
% Load magnitude is just indicative and must be tuned to correctly activate
% the non linear model.

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Mesh/Cube.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%------------------------------ ELEMENTS -----------------------------
%
GaussPts = Gauss(12,2,3);
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'Materials/Solid.dat';
%
% Create an object of the Materials class and read the materials file
mat = tabMaterials(model,topology,fileName);


%----------------------------- DOF MANAGER -----------------------------
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["BCs/bottom_fix.dat","BCs/load_z.dat","BCs/load_x.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid,dofmanager);

% Create a State object that store the solution vectors during solution
% loop
resState = State(model,grid,mat,GaussPts);

% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat','printOff');
%
% Print model initial state
printUtils.printState(resState);
%
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,linSyst,GaussPts);

% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
