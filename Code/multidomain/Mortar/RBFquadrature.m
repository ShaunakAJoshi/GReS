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
    vecPts    % pointer of each element to location in ptsRBF
  end

  properties (Access = private)
    % temporary properties for element-based sort out process
    idSlave = 0
    tempGPloc
    tempdJw
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

    function [Ns,Nm,Nmult,NbubSlave,NbubMaster] = getMortarBasisFunctions(obj,is,im)

      elemSlave = obj.mortar.getElem(2,is);

      if obj.idSlave ~= is
        % new slave element
        obj.idSlave = is;
        obj.tempdJw = getDerBasisFAndDet(elemSlave,is);
        obj.tempNs = getBasisFinGPoints(elemSlave);
        obj.tempNmult = computeMultiplierBasisF(obj.mortar,is,obj.tempNs);
        obj.tempGPloc = getGPointsLocation(elemSlave,is);
        if nargout > 3
          obj.tempNbubble = getBubbleBasisFinGPoints(elemSlave);
        end
      else
        % old slave element, sort out gp
        if all(obj.suppFlag)
          % gp already projected on all previous elements
          [Ns,Nm,Nmult,NbubSlave,NbubMaster] = deal([]);
          return
        end

        obj.tempNs = obj.tempNs(~obj.suppFlag,:);
        obj.tempNmult = obj.tempNmult(~obj.suppFlag,:);
        if nargout > 3
          obj.tempNbubble = obj.tempNbubble(~obj.suppFlag,:);
        end
        obj.tempGPloc = obj.tempGPloc(~obj.suppFlag,:);
        obj.tempdJw = obj.tempdJw(~obj.suppFlag);
      end

      % get master basis and gp in slave support
      if nargout < 5
        [Nm,obj.suppFlag] = getMasterBasisF(obj,im);
      else
        [Nm,obj.suppFlag,NbubMaster] = getMasterBasisF(obj,im);
        NbubMaster = NbubMaster(obj.suppFlag,:);
      end

      if ~any(obj.suppFlag)
        % no detected intersection
        [Ns,Nm,Nmult,NbubSlave,NbubMaster] = deal([]);
        return
      end

      Nm = Nm(obj.suppFlag,:);

      % return basis function in active gp
      Ns = obj.tempNs(obj.suppFlag,:);
      Nmult = obj.tempNmult(obj.suppFlag,:);
      if nargout > 3
        NbubSlave = obj.tempNbubble(obj.suppFlag,:);
      end
    end

    function [Nm,id,Nb] = getMasterBasisF(obj,idMaster)
      % return interpolated master basis function into slave domain
      tol = 1e-4;
      posGP = obj.tempGPloc;
      v = obj.vecPts(idMaster,:);
      nN = getElem(obj.mortar,1,idMaster).nNode;
      ptsInt = obj.ptsRBF(1:v(2),3*v(3)-[2 1 0]);
      [fiNM,id1] = obj.computeRBFfiNM(ptsInt,posGP);
      Nm = (fiNM*obj.wF(:,v(1)+1:v(1)+nN))./(fiNM*obj.w1(:,v(3)));
      Nsupp = Nm(:,[1 2 3]);
      % automatically detect supports computing interpolant
      id = all([Nsupp >= 0-tol id1],2);
      if nargout > 2
        % interpolated bubble basis function
        Nb = (fiNM*obj.wB(:,v(3)))./(fiNM*obj.w1(:,v(3)));
      end
    end

    function mat = integrate(obj,func,varargin)
      % check input
      assert(nargin(func)==numel(varargin),['Number of specified input (%i)' ...
        'not matching the integrand input (%i)'],numel(varargin),nargin(func));
      size3 = cellfun(@(x) size(x, 3), varargin);
      if ~all(size3 == size3(1))
        error('All inputs must have the same size along dimension 3.');
      end
      dJWeighed = obj.tempdJw(obj.suppFlag);
      mat = func(varargin{:});
      mat = mat.*reshape(dJWeighed,1,1,[]);
      mat = sum(mat,3);
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

      elem = obj.mortar.elements(1);
      msh = obj.mortar.mesh.msh(1);

      cells = obj.mortar.mesh.getActiveCells(1);

      obj.vecPts = zeros(max(cells),3);

      numPtsQ = (obj.nInt)^2;
      numPtsT = sum(1:obj.nInt);
      if isempty(getElement(elem,9))
        numPts = numPtsT;
      else
        numPts = numPtsQ;
      end
      
      N = sum(msh.surfaceNumVerts(cells));
     
      weighF = zeros(numPts,N);
      weigh1 = zeros(numPts,obj.mortar.mesh.nEl(1));
      weighB = weigh1;
      pts = zeros(numPts,obj.mortar.mesh.nEl(1)*3);

      k = 0;
      for i = 1:obj.mortar.mesh.nEl(1)
        im = obj.mortar.mesh.activeCells{1}(i);
        [f, ptsInt] = computeMortarBasisF(obj,im);
        nptInt = size(ptsInt,1);
        nN = getElem(obj.mortar,1,im).nNode;
        if nN < 5
          bf = computeMortarBubbleBasisF(obj,im);
        else
          bf =  ones(size(ptsInt,1),1);
          % to provisionally cope with bubble in quad9
        end
        fiMM = obj.computeRBFfiMM(ptsInt);
        % solve local system to get weight of interpolant
        warning('off','MATLAB:nearlySingularMatrix')
        x = fiMM\[f ones(size(ptsInt,1),1) bf];
        weighF(1:nptInt,k+1:k+nN) = x(:,1:nN);
        weigh1(1:nptInt,i) = x(:,end-1);
        weighB(1:nptInt,i) = x(:,end);
        pts(1:nptInt,[3*i-2 3*i-1 3*i]) = ptsInt;
        obj.vecPts(im,1) = k;
        obj.vecPts(im,2) = nptInt;
        obj.vecPts(im,3) = i;
        k = k+nN;
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
      elem = obj.mortar.getElem(1,id);
      % place interpolation points in a regular grid
      intPts = getInterpolationPoints(obj,elem);
      bf = computeBasisF(elem,intPts);
      % get coords of interpolation points in the real space
      pos = bf*coord;
    end

    function bf = computeMortarBubbleBasisF(obj,id)
      % place interpolation points in a regular grid
       elem = obj.mortar.getElem(1,id);
      intPts = getInterpolationPoints(obj,elem);
      bf = computeBubbleBasisF(elem,intPts);
    end

    function intPts = getInterpolationPoints(obj,elem)
      switch class(elem)
        case 'Triangle'
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
        case {'Quadrilateral','QuadrilateralQuadratic'}
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

