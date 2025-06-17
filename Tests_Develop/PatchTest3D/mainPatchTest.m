% simple patch test to validate hexa 27 finite elements

close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

% Set physical models 
model = ModelType("Poromechanics_FEM");

% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create the Mesh object
topology = Mesh();

% Set the mesh input file name
fileName = 'Mesh/cube.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

% Create an object of the Materials class and read the materials file
fileName = 'materialsList.dat';
mat = Materials(model,fileName);


% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(topology,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%

% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat);

% Build a structure storing variable fields at each time step
linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_Terzaghi_tetra','flagMatFile',true);


% % Write BC files programmatically with function utility 
% F = -10; % vertical force
% setTerzaghiBC('BCs',F,topology);
% 
% % Collect BC input file in a list
% fileName = ["BCs/dirFlowTop.dat","BCs/neuPorotop.dat",...
%    "BCs/dirPoroLatY.dat","BCs/dirPoroLatX.dat","BCs/dirPoroBottom.dat"];
% %
% Create an object of the "Boundaries" class 
bound = Boundaries(fileName,model,grid);

% Print model initial state
printState(printUtils,linSyst);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,linSyst);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()

