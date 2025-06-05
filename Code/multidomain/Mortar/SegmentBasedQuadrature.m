classdef SegmentBasedQuadrature < handle
  % Implement utilities to perform segment based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    mortar
  end
  
  methods
    function obj = SegmentBasedQuadrature(mortar)
      obj.mortar = mortar;
      obj.msh = obj.mortar.mesh.msh;
    end
    
    function [xiSlave,xiMaster] = segmentBasedCouple(obj,elSlave,elMaster)
      % output: nGx2xnT matrices of reference coordinate in slave and
      % master side to perform segment based integration
      % nT: numb. of triangles from delaunay triangulation on clipping
      % polygon

      % compute auxiliary plane for integration
      [P0,nP] = computeAuxiliaryPlane(obj,elSlave);
      % compute 2D coordinates of projected nodes on the aux. plane
      [coordS,coordM] = projectNodes(obj,P0,nP,elSlave,elMaster);
      % obtain intersection of polygon
      [clipX,clipY] = polyclip(coordS(:,1),coordS(:,2),coordM(:,1),coordM(:,2),1);
      assert(numel(clipX)==1,['Non unique clip polygon for master/slave pair' ...
        ' %i/%i \n'],elMaster,elSlave)
      % perform delaunay triangulation on clip polygon
      clipX = clipX{:}; clipY = clipY{:};
      triTop = delaunay(clipX,clipY);
      % project gauss points on each triangular cell into slave and master
      % side
      [xiSlave] = projectBack(obj,triTop,clipX,clipY,coordS);
      [xiMaster] = projectBack(obj,triTop,clipX,clipY,coordM);
    end
    
    function [P,n] = computeAuxiliaryPlane(obj,el)
      P = obj.msh(2).surfaceCentroid(el,:);
      elemType = obj.msh(2).surfaceVTKType(el);
      refCentroid = Elements.getElemCentroid(elemType);
      switch elemType
        
      end
    end

    function projectNodes(obj,P,n,eS,eM)
    end
  end
end

