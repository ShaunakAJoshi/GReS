function setBCfiles(leftMesh,rightMesh)
% writing BC files
if isfolder('BCs')
   [~,~,~] = rmdir('BCs','s');
   % directory become empty
end
mkdir BCs 
% custom BCs
writeBCfiles('BCs/leftLat','SurfBC','Dir','Poro',["x","y"],'left_lateral_roller',0,0,leftMesh,1); % left block lateral fix
writeBCfiles('BCs/leftBot','SurfBC','Dir','Poro',["y","z"],'left_bottom_fix',0,0,leftMesh,3); % left block bottom fix
writeBCfiles('BCs/rightBot','SurfBC','Dir','Poro',["y","z"],'right_bottom',0,0,rightMesh,4); % right block bottom fix
writeBCfiles('BCs/rightLatLoad','SurfBC','Neu','Poro',"x",'right_latLoad',[0 5 10 15 20],[0 -5 -5 0 0],rightMesh,2); % right block bottom fix
writeBCfiles('BCs/rightTopLoad','SurfBC','Neu','Poro',"z",'right_topLoad',[0 5 10 15 20],[0 0 -18 -18 0],rightMesh,3); % right block bottom fix
end

