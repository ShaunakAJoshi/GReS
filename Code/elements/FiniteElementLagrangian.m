classdef (Abstract) FiniteElementLagrangian < handle

  properties
    mesh
    indB
    indBbubble
    GaussPts
    Jref      % basis function ref. gradient
    Nref      % basis function ref.
    Jb        % bubble basis function ref. gradient
    Nb        % bubble basis function ref.
    interpOrd
    gaussOrd
    detJ
  end

  methods (Access = public)

    % Abstract class constructor
    function obj = FiniteElementLagrangian(mesh,interpolationOrder,gaussOrder)
      obj.interpOrd = interpolationOrder;
      obj.gaussOrd = gaussOrder;
      obj.mesh = mesh;
      setElement(obj);
    end

    getBasisFinGPoints(obj)
    getDerBasisFAndDet(obj)
    computeProperties(obj)
  end

  methods (Access = protected)
    setElement(obj)
    findLocBasisF(obj)
    findLocDerBasisF(obj)
  end

  methods (Static)
    function setStrainMatrix(elem)
      Nb = elem.nNode*elem.GaussPts.nNode;
      elem.indB = Poromechanics.setStrainMatIndex(Nb);
      %
      Nb = elem.nFace*elem.GaussPts.nNode;
      elem.indBbubble = Poromechanics.setStrainMatIndex(Nb);
    end
  end
end
