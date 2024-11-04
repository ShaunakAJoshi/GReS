function writeBCfiles(fName,item,type,physic,dir,bcName,time,vals,varargin)
% utility to write BC input files (only constant BC for each time step are allowed)
% Entities are identified by surface ID or by direct assignment
% Example
% writeBCfiles('SurfBC','Dir','Flow','fix_bot','1',[0 5 10],[1 2 3],msh,1)
assert(length(time)==length(vals),'Bc times and value sets must have equal size');
% writing general file
fID = fopen(strcat(fName,'.dat'),'w');
fprintf(fID,'%s            %% BC item \n',item);
if ~strcmp(item,'VolumeForce')
   fprintf(fID,'%s            %% BC type \n',type);
end
fprintf(fID,'%s            %% Physics \n',physic);
if strcmp(type,'Neu')
   assert(numel(dir)==1,'Only one direction is allowed for Neumann BCs');
   fprintf(fID,'%s \n',dir);
end
fprintf(fID,'%s            %% BC name \n',bcName);
listName = strcat(fName,'/list');
fprintf(fID,'%s \n',listName);
for i = 0:length(time)-1
   fprintf(fID,'%2.6f %s/time%i.dat \n',time(i+1),fName,i);
end
fprintf(fID,'End');

if(strcmp(physic,'Flow'))
   dir = [];
end

if ~isfolder(fName)
   mkdir(fName);
end

% writing BC list
fList = fopen(listName,'w');
if length(varargin) < 2 % direct assignment
   list = varargin{1};
else
   surf = find(ismember(varargin{1}.surfaceTag,varargin{2}));
   switch item
      case 'NodeBC'
         nodes = unique(varargin{1}.surfaces(surf,:));
         list = nodes;
      case 'SurfBC'
         list = surf;
   end
end

if isempty(dir)
   fprintf(fList,'%i         %% Number of fixed entities \n',length(list));
   fprintf(fList,'%i \n',list);
else
   if strcmp(type,'Dir')
      tmp = ismember(["x","y","z"],dir);
      fprintf(fList,'%i ',tmp*length(list));
      fprintf(fList,'   %% Number of fixed entities \n');
      list = repmat(list,sum(tmp),1);
      fprintf(fList,'%i \n',list);
   elseif strcmp(type,'Neu')
      fprintf(fList,'%i ',length(list));
      fprintf(fList,'   %% Number of fixed entities \n');
      fprintf(fList,'%i \n',list);
   end
end

% writing BC vals for each time step
for i = 1:length(time)
   t_name = strcat(fName,'/time',num2str(i-1),'.dat');
   ft = fopen(t_name,'w');
   fprintf(ft,'%%Time %2.4f \n',time(i));
   fprintf(ft,'%1.6e \n',repelem(vals(i),length(list)));
end
end
