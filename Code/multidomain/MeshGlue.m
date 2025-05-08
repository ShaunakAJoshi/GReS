classdef MeshGlue < handle
  %
  properties
    idMaster
    idSlave
    dofMaster % dof id of master entities (wrt local domain)
    dofSlave % dof id of slave entities (wrt local domain)
    slaveMat % cross grid matrix mapping slave dofs to multipliers
    masterMat % cross grid matrix mapping master dofs to multipliers
  end
  
  methods
    function obj = MeshGlue(inputStruct,domains)
      % 
      
    end
    
    function outputArg = method1(obj,inputArg)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      outputArg = obj.Property1 + inputArg;
    end
  end

  methods (Static)
    function interfaceStruct = buildInterfaceStruct(fileName,modelStruct)
      % read interface file and construct array of MeshGlue objects
      interfStr = readstruct(fileName);
      interfStr = interfStr.Interface;
      nInterfaces = numel(interfStr);
      interfaceStruct(nInterfaces, 1) = MeshGlue();  % Preallocate array
      for i = 1:nInterfaces
        idMaster = interfStr(i).Master.idAttribute;
        idSlave = interfStr(i).Slave.idAttribute;
        type = interfStr(i).Type;
        if strcmp(type,'MeshTying')
          interfaceStruct(i) = MeshGlue(interfStr(i), ...
            modelStruct([idMaster,idSlave]));
        elseif strcmp(type,'Fault')
          % not yet implemented!
          interfaceStruct(i) = Fault(interfStr(i), ...
            modelStruct([idMaster,idSlave]));
        else
          error(['Invalid interface law type for interface %i in file' ...
            '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
      end
    end
  end
end

