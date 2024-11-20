clear
close all

% Get the full path of the script
scriptFullPath = mfilename('fullpath'); 
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);


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
fileName = 'Mesh/Block.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'Materials/MaterialsList.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);
%
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
dofmanager = DoFManager(topology,model);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
%
% apply BCs
writeBCfiles('BCs/Bot','SurfBC','Dir','Poro',["y","z"],'bottom',0,0,topology,3); % left block lateral fix
writeBCfiles('BCs/LatBound','SurfBC','Dir','Poro',["x","y"],'left_bottom_fix',0,0,topology,1); % left block bottom fix
writeBCfiles('BCs/LatLoad','SurfBC','Neu','Poro',"x",'right_latLoad',[0 5 10 15 20],[0 -5 -5 0 0],topology,2); % right block bottom fix
%
% Set the input file
fileName = ["BCs/Bot.dat","BCs/LatBound.dat","BCs/LatLoad.dat"];
bound = Boundaries(fileName,model,grid,dofmanager);


resState = State(model,grid,mat,GaussPts);

% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat','printOff');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,linSyst,GaussPts);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()



