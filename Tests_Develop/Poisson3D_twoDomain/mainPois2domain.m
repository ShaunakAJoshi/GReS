clear
close all

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));

%% model commons
simParam = SimulationParameters('simParam.dat');
% base structure to write xml file
strDomain = readstruct('Domains/domain2block.xml');
interfFile = 'interface.xml';

%% INPUT
% base domain size
N_0_l = 4;
N_0_r = 5;
% number of refinement
nref = 1;
[h,L2,H1] = deal(zeros(nref,1));

%% convergence loop
for i = 1:nref
  N_i_l = N_0_l*2^(i-1);
  N_i_r = N_0_r*2^(i-1);

  fprintf('Running mesh refinement %i \n',i);

  % run script to get refined mesh
  fname = strcat('domain_',num2str(i));
  command = "python Mesh/domain.py "  + fname...
    + " " + num2str(N_i_l) + " " + num2str(N_i_r);
  system(command);
  
  meshFile = char(fname+".vtk");
  mesh = Mesh();
  mesh.importMesh(meshFile);

  % set up bc file
  nodes = unique(mesh.surfaces(ismember(mesh.surfaceTag,[2 4]),:));
  c = mesh.coordinates(nodes,:);
  vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
  vals = reshape(vals,[],1);
  writeBCfiles('bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);

  clear mesh

  % write mesh to domain file
  strDomain.Domain.Geometry = fullfile(fname+".vtk");
  domainFile = fullfile('Domains','domain2block.xml');
  writestruct(strDomain,domainFile);
 
  % processing Poisson problem
  domains = buildModelStruct_new(domainFile,simParam);
  domains.Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);
  
  [interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);
  % set up analytical solution
  
  solver = MultidomainFCSolver(simParam,domains,interfaces);
  solver.NonLinearLoop();
  solver.finalizeOutput();

  %runPoisson;

  pois = getSolver(domains.Discretizer,'Poisson');
  [L2(i),H1(i)] = pois.computeError_v2();
  h(i) = 1/N_i_r;
  fprintf('Max absolute error is: %1.6e \n',max(abs(pois.state.data.err)));
end

% compute convergence order
L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));

%% plotting convergence profiles
% figure(1)
% loglog(h,L2,'-ro')
% hold on
% loglog(h,H1,'-b^')
% xlabel('h')
% ylabel('error_norm')
% legend('L2','H1')




