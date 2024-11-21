function  plotFunction(mesh, foldName, funct, name, varargin)
% shortcut function to produce VTK output
% array solutions can be handled

outVTK = VTKOutput(mesh, foldName); % create VTK object
ncomp = size(funct,2);
assert(numel(name)==ncomp,'Wrong array size for input name');

pointData = repmat(struct('name',[],'data',[]),ncomp,1);
for i = 1:size(funct,2)
    pointData(i).name = convertStringsToChars(name(i));
    pointData(i).data = funct(:,i);
end
if isempty(varargin) || strcmpi(varargin{1},'node')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(0, pointData, [], [], []);
   else
      outVTK.writeVTKFile(0, [], [], pointData, []);
   end
elseif strcmpi(varargin{1},'elem')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(0, [], pointData, [], []);
   else
      outVTK.writeVTKFile(0, [], [], [], pointData);
   end
end
   outVTK.finalize()
end

