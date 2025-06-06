classdef Triangle < FiniteElementLagrangian
  % TRIANGLE element class

  properties (Constant)
    centroid = [0.33,0.33]
    coordLoc = [0 0; 1 0; 0 1];
    vtkType = 5
    nNode = 3
    nFace = 1
  end



  methods (Access = public)

    function [mat] = getDerBasisF(obj,el)
      % compute derivatives of the basis functions for element in real
      % space
      inv_A = inv([1 obj.mesh.coordinates(obj.mesh.surfaces(el,1),1:2);
        1 obj.mesh.coordinates(obj.mesh.surfaces(el,2),1:2);
        1 obj.mesh.coordinates(obj.mesh.surfaces(el,3),1:2)]);
      mat = inv_A(2:3,:);
    end

    function dN = computeDerBasisF(obj,varargin)
      dN = obj.J1;
    end

     function [outVar1,outVar2] = getDerBasisFAndDet(obj,in)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % way to call this method: if el is a scalar (element idx) the 3D
      % coordinates are retrieved by the corresponding mesh object. Only
      % the determinant is returned
      % if in is not scalar, it is a 4x2 list of 2 coordinates, the
      % gradient matrix and the determinant are returned

      if isscalar(in)
        % 3D setting
        coord = obj.mesh.coordinates(obj.mesh.surfaces(in,:),1:obj.mesh.nDim);
        J = pagemtimes(obj.Jref,coord);
        for i = 1:obj.GaussPts.nNode
          obj.detJ(i) = norm(cross(J(1,:,i),J(2,:,i)),2);
        end
        outVar1 = obj.detJ.*(obj.GaussPts.weight)';
      else
        % 2D setting: in is a given list of x-y coordinates
        J = pagemtimes(obj.Jref,in);
        for i=1:obj.GaussPts.nNode
          J(:,:,i) = inv(J(:,:,i));
          obj.detJ(i) = det(J(:,:,i));
        end
        outVar1 = pagemtimes(J,obj.Jref);
        outVar2 = obj.detJ.*(obj.GaussPts.weight)';
      end
    end

    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.Nref;
    end

    function N = computeBasisF(obj, list)
      % Find the value the basis functions take at some  reference points defined in
      % a list
      if size(list,1) > 1
        N = zeros(size(list,1),3);
        N(:,1) = 1-list(:,1)-list(:,2);
        N(:,2) = list(:,1);
        N(:,3) = list(:,2);
      else
        N = zeros(1,3);
        N(1) = 1-list(1)-list(2);
        N(2) = list(1);
        N(3) = list(2);
      end
    end

    function gPCoordinates = getGPointsLocation(obj,el)
      % Get the location of the Gauss points in the element in the physical
      % space
      gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
    end

    function [area,cellCentroid] = findAreaAndCentroid(obj,idTri)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idTri),1);
      cellCentroid = zeros(length(idTri),3);
      i = 0;
      for el = idTri'
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        coord = obj.mesh.coordinates(obj.mesh.surfaces(idTri,:),:);
        cellCentroid(i,:) = 1/3*(sum(coord,1));
      end
    end

    function n = computeNormal(obj,idTri)
      % compute normal vector of cell idQuad
      n = zeros(length(idTri),3);
      for el = idTri
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
        a = findAreaAndCentroid(obj,i);
        n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + a/3;
      end
      n_a = n_a(unique(surfMsh.surfaces));
    end


    function computeProperties(obj)
      % update the parent mesh object computing cell properties
      idTri = find(obj.mesh.surfaceVTKType == obj.vtkType);
      [area,centr] = findAreaAndCentroid(obj,idTri);
      obj.mesh.surfaceCentroid(idTri,:) = centr;
      obj.mesh.surfaceArea(idTri,:) = area;
    end
  end

  methods (Access = protected)

    function findLocBasisF(obj)
      % Find the value the basis functions take at the Gauss points
      Ntmp = zeros(obj.GaussPts.nNode,3);
      Ntmp(:,1) = 1-obj.GaussPts.coord(:,1)-obj.GaussPts.coord(:,2);
      Ntmp(:,2) = obj.GaussPts.coord(:,1);
      Ntmp(:,3) = obj.GaussPts.coord(:,2);
      obj.Nref = Ntmp;
      if obj.GaussPts.nNode == 1
        obj.Nref = obj.Nref';
      end
    end

    function findLocDerBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      obj.Jref = [-1 1 0; -1 0 1; 0 0 0];
    end

    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.gaussOrd);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
    end

  end
end
