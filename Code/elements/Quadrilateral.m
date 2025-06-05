 classdef Quadrilateral < handle
  % QUADRILATERAL element class

  properties (Access = public)
%     vol
%     volNod
%     cellCentroid
  end

  properties (Constant)
    centroid = [0,0]
    coordLoc = [-1 -1;
      1 -1;
      1  1;
      -1  1;]

  end

  properties (Access = private)   % PRIVATE
%
% NODE ORDERING ASSUMPTION (same as Gmsh output):
% Quadrilatersl (4 nodes):             
% 
%       | v
% 4-----|-----3            
% |     |     |           
% |     |     |           
% |     +-----|---->                
% |           |   u  
% |           |      
% 1-----------2    


    
    GaussPts
    J1
    mesh
    J
    detJ
    N1
    Nb
    Jb
  end

  methods (Access = public)

    % Class constructor method
    function obj = Quadrilateral(msh,GPoints)
       obj.setQuad(msh,GPoints);
    end

      

    function [outVar1,outVar2] = getDerBasisFAndDet(obj,in)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and detJ
      % way to call this method: if el is a scalar (element idx) the 3D
      % coordinates are retrieved by the corresponding mesh object. Only
      % the determinant is returned
      % if in is not scalar, it is a 4x2 list of 2 coordinates, the
      % gradient matrix and the determinant are returned

      if isscalar(in)
        % 3D setting
        coord = obj.mesh.coordinates(obj.mesh.surfaces(in,:),1:obj.mesh.nDim);
        obj.J = pagemtimes(obj.J1,coord);
        for i = 1:obj.GaussPts.nNode
          obj.detJ(i) = norm(cross(obj.J(1,:,i),obj.J(2,:,i)),2);
          outVar1 = obj.detJ.*(obj.GaussPts.weight)';
        end
      else
        % 2D setting: in is a given list of x-y coordinates
        obj.J = pagemtimes(obj.J1,in);
        for i=1:obj.GaussPts.nNode
          obj.J(:,:,i) = inv(obj.J(:,:,i));
          obj.detJ(i) = det(obj.J(:,:,i));
        end
        outVar1 = pagemtimes(obj.J,obj.J1);
        outVar2 = obj.detJ.*(obj.GaussPts.weight)';
      end
    end

    
    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.N1;
    end

    function NbMat = getBubbleBasisFinGPoints(obj)
      NbMat = obj.Nb;
      NbMat = reshape(NbMat,obj.GaussPts.nNode,[]);
    end

 

    function [area,cellCentroid] = findAreaAndCentroid(obj,idHexa)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idHexa),1);
      cellCentroid = zeros(length(idHexa),3);
      i = 0;
      for el = idHexa
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el,3);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        gPCoordinates = getGPointsLocation(obj,el);
        cellCentroid(i,:) = obj.detJ * gPCoordinates/area(i);
      end
    end
    


    function nodeArea = findNodeArea(obj,idQuad)
      nodeArea = zeros(4*length(idQuad),1);
      ptr = 0;
      for el = idQuad
        dJWeighed = obj.getDerBasisFAndDet(el,3);
        nodeArea(ptr+1:ptr+4) = obj.N1'*dJWeighed';
        ptr = ptr + 4;
      end
    end



    function n = computeNormal(obj,idQuad)
        % compute normal vector of cell idQuad
        n = zeros(length(idQuad),3);
        for el = idQuad
            nodeCoord = obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
            v1 = nodeCoord(1,:) - nodeCoord(2,:);
            v2 = nodeCoord(2,:) - nodeCoord(3,:);
            n(el,:) = cross(v1,v2);
            n(el,:) = n(el,:)/norm(n(el,:));
        end
    end



    function n_a = computeAreaNod(obj,surfMsh)
       % compute area associated to each node of a surface mesh
       n_a = zeros(max(surfMsh.surfaces,[],'all'),1);
       for i = 1:length(surfMsh.surfaces)
          n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + findNodeArea(obj,i);
       end
       n_a = n_a(unique(surfMsh.surfaces));
    end



    function gPCoordinates = getGPointsLocation(obj,el)
        % Get the location of the Gauss points in the element in the physical
        % space
        gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
    end



    function N = computeBasisF(obj, coord)
      % Find the value the basis functions take at some  reference points
      % whose 2D coordinates are store in coord
      N = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*coord(i,1)).* ...
        (1+obj.coordLoc(j,2).*coord(i,2)), ...
        (1:size(coord,1))',1:obj.mesh.surfaceNumVerts(1));
      if size(N,2) ~= obj.mesh.surfaceNumVerts(1)
        N = N';
      end
    end

    function N = computeBubbleBasisF(obj, coord)
      % Find the value the bubble basis functions take at some  reference
      % points whose 2D coordinates are store in coord
      N = arrayfun(@(i) (1-coord(i,1)^2).*(1-coord(i,2)^2),(1:size(coord,1)));
      N = N';
    end



    function dN = computeDerBasisF(obj, list)
        % Compute derivatives in the reference space for all Gauss points
        % d(N)/d\csi
        d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
            (1+obj.coordLoc(j,2).*list(i,2)), ...
            (1:size(list,1)),1:obj.mesh.surfaceNumVerts(1));
        %
        % d(N)/d\eta
        d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
            (1+obj.coordLoc(j,1).*list(i,1)), ...
            (1:size(list,1)),1:obj.mesh.surfaceNumVerts(1));
        %
        dN = [d1';d2'];
    end

  end

  methods (Access = private)

    function findLocDerBasisF(obj,varargin)
      % Compute derivatives in the reference space for all Gauss points
      obj.J1 = zeros(2,obj.mesh.surfaceNumVerts(1),obj.GaussPts.nNode);
      %
      % d(N)/d\csi
      if isempty(varargin)
        d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
        %
        % d(N)/d\eta
        d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
        % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
        %
        obj.J1(1,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d1';
        obj.J1(2,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d2';
      else
        refCoord = varargin{1};
        d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
          (1+obj.coordLoc(j,2).*refCoord(i,2)), ...
          (1:size(refCoord,1))',1:obj.mesh.surfaceNumVerts(1));
        %
        % d(N)/d\eta
        d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
          (1+obj.coordLoc(j,1).*refCoord(i,1)), ...
          (1:size(refCoord,1))',1:obj.mesh.surfaceNumVerts(1));
        % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
        %
        obj.J1(1,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d1';
        obj.J1(2,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d2';
      end
    end
    


    function findLocBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      if isempty(varargin)
        obj.N1 = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
        if obj.GaussPts.nNode == 1
          obj.N1 = obj.N1';
        end
      else
        % compute basis at given reference point (xi,eta)
        refCoord = varargin{1};
        bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*refCoord(i,1)).* ...
          (1+obj.coordLoc(j,2).*refCoord(i,2)), ...
          (1:size(refCoord,1))',1:obj.mesh.surfaceNumVerts(1));
      end
    end


    function findLocBubbleBasisF(obj)
      % Find the value the basis functions take at the Gauss points

      g = obj.GaussPts.coord;
      bub = @(x,y) 1-g(x,y)^2;

      obj.Nb = zeros(obj.GaussPts.nNode,1);
      obj.Nb =  arrayfun(@(i) bub(i,1).*bub(i,2),1:obj.GaussPts.nNode);
    end    


    function setQuad(obj,msh,GPoints)
      obj.mesh = msh;
      obj.GaussPts = GPoints;
      findLocDerBasisF(obj);
      findLocBasisF(obj);
      findLocBubbleBasisF(obj);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
    end
  end

  end
