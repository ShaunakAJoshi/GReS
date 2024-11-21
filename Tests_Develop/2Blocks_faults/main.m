clear
close all

% Get the full path of the script
scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('2blocks model \n')
fprintf('___________________\n\n')
[leftMesh,rightMesh] = deal(Mesh(),Mesh());

simulTag = 'tetra';
switch simulTag
   case 'tetra'
      domainFile = 'domains_tetra.dat';
      meshLeftFile = 'Mesh/LeftBlock_tetra.msh';
      meshRightFile = 'Mesh/RightBlock_tetra.msh';
   case 'hexa'
      domainFile = 'domains_hexa.dat';
      meshLeftFile = 'Mesh/LeftBlock_hexa.msh';
      meshRightFile = 'Mesh/RightBlock_hexa.msh';
end

leftMesh.importGMSHmesh(meshLeftFile);
rightMesh.importGMSHmesh(meshRightFile);

% write BC files
setBCfiles(leftMesh,rightMesh);

plotFunction(leftMesh,'out_L',zeros(leftMesh.nNodes,1),"sol");
plotFunction(rightMesh,'out_R',zeros(rightMesh.nNodes,1),"sol");


simParam = SimulationParameters('simParam.dat');
model = buildModelStruct(domainFile,simParam);
mG = MeshGlue(model,'interface.dat');

% validate mortar operator construction
f = @(y,z) sin(y)+cos(z);
cs = mG.interfaces.mortar.intSlave.coordinates;
cm = mG.interfaces.mortar.intMaster.coordinates;
fSAnal = f(cs(:,2),cs(:,3));
fMaster = f(cm(:,2),cm(:,3));
fSInterp = mG.interfaces.InterpOperator*fMaster;
e = fSInterp-fSAnal;
e = sqrt(sum(e.^2));

% fault material parameters
coes = 0;
phi = 30; % degrees

solver = mortarFaults(simParam,mG,coes,phi);

