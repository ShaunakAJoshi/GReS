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
    nGtri
    elems     % provisional elems instances for each master/slave pair
  end
  
  methods
    function obj = SegmentBasedQuadrature(mortar,nGtri)
      obj.mortar = mortar;
      obj.msh = obj.mortar.mesh.msh;
      obj.nGtri = nGtri;
    end

    function [Ns,Nm,Nmult,Nbub] = getMortarBasisFunctions(obj,is,im)
      [xiSlave,xiMaster,dJwTri] = segmentBasedCouple(obj,is,im);
      nTri = size(dJwTri,2);
      elemSlave = obj.mortar.getElem(2,is);
      elemMaster = obj.mortar.getElem(1,im);
      Nm = zeros(obj.nGtri, elemMaster.nNode,nTri);
      Ns = deal(zeros(obj.nGtri,elemSlave.nNode,nTri));
      Nbub = deal(zeros(obj.nGtri,1,nTri));
      Nmult = ones(obj.nGtri,1,nTri);               % for P0 elements!
      for i = 1:nTri
        xiM = xiMaster(:,:,i);
        xiS = xiSlave(:,:,i);
        Nm(:,:,i) = elemMaster.computeBasisF(xiM);
        Nbub(:,:,i) = elemSlave.computeBubbleBasisF(xiS);
        Ns(:,:,i) = elemSlave.computeBasisF(xiS);
      end
    end
    
    function [xiSlave,xiMaster,dJwTri] = segmentBasedCouple(obj,elSlave,elMaster)
      % output: nGx2xnT matrices of reference coordinate in slave and
      % master side to perform segment based integration
      % dJwTri: weighed jacobian determinant for each pallet. size nGxnTri
      % nT: numb. of triangles from delaunay triangulation on clipping
      % polygon
      % compute auxiliary plane for integration
      obj.elems = [getElem(obj.mortar,1,elMaster),...
                   getElem(obj.mortar,2,elSlave)];
      [P0,nP] = computeAuxiliaryPlane(obj,elSlave);
      % compute 2D coordinates of projected nodes on the aux. plane
      [coordS,coordM] = projectNodes(obj,P0,nP,elSlave,elMaster);
      % obtain intersection of polygon
      [clipX,clipY] = polyclip(coordS(:,1),coordS(:,2),coordM(:,1),coordM(:,2),1);
      assert(numel(clipX)==1,['Non unique clip polygon for master/slave pair' ...
        ' %i/%i \n'],elMaster,elSlave)
      % perform delaunay triangulation on clip polygon
      % assumption: only one clip polygon results from intersection
      coordClip = [clipX{:} clipY{:}];
      topolClip = delaunay(coordClip(:,1),coordClip(:,2));
      % project gauss points on each triangular cell into slave and master
      % side
      [xiSlave] = projectBack(obj,2,topolClip,coordClip,coordS);
      [xiMaster] = projectBack(obj,1,topolClip,coordClip,coordM);
      tri = Triangle(1,obj.nGtri);
      nTri = size(topolClip,1);
      dJwTri = zeros(obj.nGtri,nTri);
      for i = 1:nTri
        triVert = coordClip(topolClip(i,:),:);
        dJwTri(:,i) = getDerBasisFAndDet(tri,triVert);
      end
    end
    
    function [P,n] = computeAuxiliaryPlane(obj,el)
      P = obj.msh(2).surfaceCentroid(el,:);
      n = obj.elems(2).computeNormal(el,obj.elems(2).centroid);
    end

    function xi = projectBack(obj,side,topolTri,clipCoord,elemCoord)
      % return reference coordinates in master/slave space for GP in
      % triangle facets after intersection
      % output is a 3D matrices of size nGx2xnTri
      nTri = size(topolTri,1);
      xi = zeros(obj.nGtri,2,nTri);
      % netwon params
      itMax = 10;
      tol = 1e-9;
      tri = Triangle(1,obj.nGtri);                % define reference triangle
      for i = 1:nTri
        coordTri = clipCoord(topolTri(i,:),:);
        coordGPtri = getGPointsLocation(tri,coordTri);
        for g = 1:obj.nGtri
          rhs = (obj.elems(side).computeBasisF(xi(g,:))*elemCoord)' - coordGPtri(g,:)';
          iter = 0;
          while (norm(rhs,2) > tol) && (iter < itMax)
            iter = iter+1;
            J = (obj.elems(side).computeDerBasisF(xi(g,:,i))*elemCoord)';
            dxi = J\(-rhs);
            xi(g,:,i) = xi(g,:,i) + dxi';
            rhs = (obj.elems(side).computeBasisF(xi(g,:,i))*elemCoord)' - coordGPtri(g,:)';
          end
          fl = FiniteElementLagrangian.checkInRange(obj.elems(side),xi(g,:,i));
          assert(fl,['GP %i in triangle %i out of range of coordinates of' ...
            'reference space']);
        end
      end
    end

    function [xS,xM] = projectNodes(obj,P,n,eS,eM)

      % get 2D direction of plane
      % Choose arbitrary vector not parallel to n
      if abs(n(1)) < 0.9
        temp = [1; 0; 0];
      else
        temp = [0; 1; 0];
      end
      % First direction in the plane
      d1 = cross(n, temp);
      d1 = d1 / norm(d1);
      % Second direction in the plane
      d2 = cross(n, d1);
      
      x = cell(2,1);
      e = [eM,eS];
      for i = 1:2
        c = obj.msh(i).coordinates(obj.msh(i).surfaces(e(i),:),:);
        cn = (c - P)*n;
        projC = c - repmat(n',obj.elems(i).nNode,1).*cn;
        % project in 2D
        x{i} = (projC - P)*[d1 d2];
      end
      [xM,xS] = deal(x{:});
      
    end
  end
end

