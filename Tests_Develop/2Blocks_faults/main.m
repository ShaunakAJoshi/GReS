clear
close all

%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('2blocks model \n')
fprintf('___________________\n\n')
[leftMesh,rightMesh] = deal(Mesh(),Mesh());

leftMesh.importGMSHmesh('Mesh/LeftBlock.msh');
rightMesh.importGMSHmesh('Mesh/RightBlock.msh');


plotFunction(rightMesh,'out',zeros(rightMesh.nNodes,1));
% writing BC files
% writeBCfiles('BCs/leftLat','SurfBC','Dir','Poro',["x","y"],'left_lateral_roller',0,0,leftMesh,1); % left block lateral fix
% writeBCfiles('BCs/leftBot','SurfBC','Dir','Poro',["y","z"],'left_bottom_fix',0,0,leftMesh,3); % left block bottom fix
% writeBCfiles('BCs/rightBot','SurfBC','Dir','Poro',["y","z"],'right_bottom',0,0,rightMesh,4); % right block bottom fix
% writeBCfiles('BCs/rightLatLoad','SurfBC','Neu','Poro',"x",'right_latLoad',[0 4 10 15 20],[0 -18 -18 0 0],rightMesh,3); % right block bottom fix
% writeBCfiles('BCs/rightTopLoad','SurfBC','Neu','Poro',"z",'right_topLoad',[5 10 15 20],[0 -5 -5 0],rightMesh,2); % right block bottom fix

simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains.dat',simParam);
mG = MeshGlue(model,'interface.dat');

