classdef RBFquadrature < handle
  % Radial Basis Function class
  % Used to interpolate master basis functions into the slave side

  properties
    mortar    % instance of mortar object (store mesh info)
    nInt      % number of interpolation points per master element
    wF        % weight of RBF for primary variable basis function
    w1        % weight of rescaling RBF
    wB        % weight of RBF for bubble basis functions
    ptsRBF    % position of interpolation points
  end

  properties
    % temporary properties for element-based sort out process
    idSlave = 0
    tempGPloc
    tempNs
    tempNmult
    tempNbubble
    suppFlag       
  end

  methods
    function obj = RBFquadrature(mortar,nInt)
      %
      obj.mortar = mortar;
      obj.nInt = nInt;
      getWeights(obj); % interpolation weights for master basis function
    end

  end
  
  methods (Access = public)

    function [Ns,Nm,Nmult,Nbub] = getMortarBasisFunctions(obj,is,im)
      elemSlave = obj.mortar.getElem(2,is);
      if obj.idSlave ~= is
        % new slave element
        obj.idSlave = is;
        obj.tempNs = getBasisFinGPoints(elemSlave);
        obj.tempNmult = ones(size(obj.tempNs,1),1);
        obj.tempNbubble = getBubbleBasisFinGPoints(elemSlave);
        obj.tempGPloc = getGPointsLocation(elemSlave,is);
      else
        % old slave element, sort out gp
        obj.tempNs = obj.tempNs(~obj.suppFlag);
        obj.tempNmult = obj.tempNmult(~obj.suppFlag);
        obj.tempNbubble = obj.tempNbubble(~obj.suppFlag);
        obj.tempGPloc = obj.tempGPloc(~obj.suppFlag);
      end
      % get master basis and gp in slave support
      [Nm,obj.suppFlag] = getMasterBasisF(obj,im);
      Nm = Nm(obj.suppFlag);
      % return basis function in active gp
      Ns = obj.tempNs(obj.suppFlag);
      Nmult = obj.tempNmult(obj.suppFlag);
      Nbub = obj.tempNbubble(obj.suppFlag);
    end

    function [Nm,id,Nb] = getMasterBasisF(obj,idMaster)
      % return interpolated master basis function into slave domain
      tol = 1e-4;
      posGP = obj.tempGPloc;
      ptsInt = obj.ptsRBF(:,repNum(3,idMaster));
      [fiNM,id1] = obj.computeRBFfiNM(ptsInt,posGP);
      Nm = (fiNM*obj.wF(:,repNum(obj.mortar.mesh.nN(1),idMaster)))./(fiNM*obj.w1(:,idMaster));
      Nsupp = Nm(:,[1 2 3]);
      % automatically detect supports computing interpolant
      id = all([Nsupp >= 0-tol id1],2);
      if nargout > 2
        % interpolated bubble basis function
        Nb = (fiNM*obj.wB(:,idMaster))./(fiNM*obj.w1(:,idMaster));
      end
    end

%     function [Nm,id] = getMasterBubbleBasisF(obj,idMaster,posGP)
%       % return interpolated master basis function into slave domain
%       tol = 1e-4;
%       ptsInt = obj.ptsRBF(:,repNum(3,idMaster));
%       [fiNM,id1] = obj.computeRBFfiNM(ptsInt,posGP);
%       Nm = (fiNM*obj.wF(:,repNum(obj.mortar.mesh.nN(1),idMaster)))./(fiNM*obj.w1(:,idMaster));
%       Nsupp = Nm(:,[1 2 3]);
%       % automatically detect supports computing interpolant
%       id = all([Nsupp >= 0-tol id1],2);
%     end
  end

  methods (Access = private)
    %
    function getWeights(obj)
      switch obj.mortar.mesh.cellType(1)
        case 10
          numPts = sum(1:obj.nInt);
        case 12
          numPts = (obj.nInt)^2;
      end
      weighF = zeros(numPts,obj.mortar.mesh.nEl(1)*obj.mortar.mesh.nN(1));
      weigh1 = zeros(numPts,obj.mortar.mesh.nEl(1));
      weighB = weigh1;
      pts = zeros(numPts,obj.mortar.mesh.nEl(1)*3);
      for i = 1:obj.mortar.mesh.nEl(1)
        [f, ptsInt] = computeMortarBasisF(obj,i);
        bf = computeMortarBubbleBasisF(obj);
        fiMM = obj.computeRBFfiMM(ptsInt);
        % solve local system to get weight of interpolant
        warning('off','MATLAB:nearlySingularMatrix')
        weighF(:,repNum(obj.mortar.mesh.nN(1),i)) = fiMM\f;
        weigh1(:,i) = fiMM\ones(size(ptsInt,1),1);
        weighB(:,i) = fiMM\bf;
        pts(:,repNum(3,i)) = ptsInt;
      end
      obj.wF = weighF;
      obj.w1 = weigh1;
      obj.wB = weighB;
      obj.ptsRBF = pts;
      % loop trough master elements and interpolate Basis function
    end

    function fiMM = computeRBFfiMM(obj,ptsInt)
      r = obj.computeRBFradius(ptsInt);
      d = sqrt((ptsInt(:,1) - ptsInt(:,1)').^2 + (ptsInt(:,2) - ptsInt(:,2)').^2 + (ptsInt(:,3) - ptsInt(:,3)').^2);
      fiMM = obj.rbfInterp(d,r);
    end

    function [fiNM,id] = computeRBFfiNM(obj,ptsInt,ptsGauss)
      % id: id of points that has a distance < r with at least one
      % master points
      d = sqrt((ptsGauss(:,1) - ptsInt(:,1)').^2 + (ptsGauss(:,2) - ptsInt(:,2)').^2 + (ptsGauss(:,3) - ptsInt(:,3)').^2);
      r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
      id = ~all(d>=r,2);
      fiNM = obj.rbfInterp(d,r);
    end

    function [bf,pos] = computeMortarBasisF(obj,id)
      % evaluate shape function in the real space and return position of
      % integration points in the real space
      surfNodes = obj.mortar.mesh.msh(1).surfaces(id,:);
      coord = obj.mortar.mesh.msh(1).coordinates(surfNodes,:);
      elem = obj.mortar.getElem(1);
      % place interpolation points in a regular grid
      intPts = getInterpolationPoints(obj);
      bf = computeBasisF(elem,intPts);
      % get coords of interpolation points in the real space
      pos = bf*coord;
    end

    function bf = computeMortarBubbleBasisF(obj)
      % place interpolation points in a regular grid
      intPts = getInterpolationPoints(obj);
      elem = obj.mortar.getElem(1);
      bf = computeBubbleBasisF(elem,intPts);
    end

    function intPts = getInterpolationPoints(obj)
      switch obj.mortar.mesh.cellType(1)
        case 10
          % uniform grid in triangle
          intPts = zeros(sum(1:obj.nInt),2);
          p = linspace(0,1,obj.nInt);
          k = obj.nInt;
          c = 0;
          for i = 1:obj.nInt
            intPts(c+1:c+k,1) = p(1:k)';
            c = c + k;
            k = k-1;
          end
          intPts(:,2) = repelem(p,obj.nInt:-1:1);
        case 12
          intPts = linspace(-1,1, obj.nInt);
          [y, x] = meshgrid(intPts, intPts);
          intPts = [x(:), y(:)];
      end
    end
  end

  methods (Static)

    function rbf = rbfInterp(d,r)
      % compute row of rbf interpolation matrix
      rbf = exp(-d.^2/r^2);
    end

    function r = computeRBFradius(ptsInt)
      % compute the radius of the RBF interpolation based on
      % coordinates of interpolation points
      r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + ...
        (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + ...
        (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
    end
  end

end

