classdef interfaceMesh
  % Manager for topological information and operation on interfaces between
  % meshes

  
  properties
    Property1
  end
  
  methods
    function obj = interfaceMesh(inputArg1,inputArg2)
      %INTERFACEMESH Construct an instance of this class
      %   Detailed explanation goes here
      obj.Property1 = inputArg1 + inputArg2;
    end
    
    function outputArg = method1(obj,inputArg)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      outputArg = obj.Property1 + inputArg;
    end
  end
end

