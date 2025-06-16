classdef ElementBasedQuadrature < handle
  % Implement utilities to perform element based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    mortar
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
    function obj = ElementBasedQuadrature(mortar)
      obj.mortar = mortar;
      obj.msh = obj.mortar.mesh.msh;
    end

    function [Ns,Nm,Nmult,NbubSlave,NbubMaster] = getMortarBasisFunctions(obj,is,im)

      elemSlave = obj.mortar.getElem(2,is);
      elemMaster = obj.mortar.getElem(1,im);

      % new slave element
      if obj.idSlave ~= is
        obj.idSlave = is;
        obj.tempdJw = getDerBasisFAndDet(elemSlave,is);
        obj.tempNs = getBasisFinGPoints(elemSlave);
        obj.tempNmult = computeMultiplierBasisF(obj.mortar,is,obj.tempNs);
        obj.tempGPloc = elemSlave.GaussPts.coord;
        obj.suppFlag = false(size(obj.tempNs,1),1);
        if nargout > 3
          obj.tempNbubble = getBubbleBasisFinGPoints(elemSlave);
        end
      else         % old slave element, sort out gp
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
      [xiProj,obj.suppFlag] = projectGP(obj,is,im);
      Nm = elemMaster.computeBasisF(xiProj);

      if nargout > 4
        NbubMaster = elemMaster.computeBubbleBasisF(xiProj);  
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

    function [xiM,id] = projectGP(obj,is,im)
      % xi: reference coordinates of the gauss point
      % get nodal normal
      supp = ~obj.suppFlag;
      elM = getElem(obj.mortar,1,im);
      elS = getElem(obj.mortar,2,is);
      nodeS = obj.msh(2).surfaces(is,:);
      coordS = obj.msh(2).coordinates(nodeS,:);
      X = elS.Nref(supp,:)*coordS;                    % real position of gauss pts
      xiS = obj.tempGPloc;
      ngp = size(xiS,1);
      xiM = zeros(ngp,2);
      id = false(ngp,1);
      itMax = 8;
      tol = 1e-9;
      for i = 1:ngp
        Ns = elS.computeBasisF(xiS(i,:));
        ng = Ns*obj.mortar.mesh.avgNodNormal{2}(nodeS,:); % slave normal at GP
        xiM(i,:) = elS.centroid;                                 % initial guess
        iter = 0;
        w = 0;
        nodeM = obj.msh(1).surfaces(im,:);
        coordM = obj.msh(1).coordinates(nodeM,:);
        Nm = elM.computeBasisF(xiM(i,:));
        rhs = (Nm*coordM - w*ng - X(i,:))';
        %
        while (norm(rhs,2) > tol) && (iter < itMax)
          iter = iter+1;
          J1 = (elM.computeDerBasisF(xiM(i,:))*coordM)';
          J = [J1 -ng'];
          ds = J\(-rhs);
          xiM(i,:) = xiM(i,:) + (ds(1:2))';
          w = w + ds(3);
          Nm =  elM.computeBasisF(xiM(i,:));
          rhs = (Nm*coordM - w*ng - X(i,:))';
        end

        % check if gp lies in master reference space
        id(i) = FiniteElementLagrangian.checkInRange(elM,xiM(i,:));
      end
    end


    function mat = integrate(obj,func,varargin)
      % element based integration
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
    
  end
end

