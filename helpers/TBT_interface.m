function [initMat, iterMat] = TBT_interface(p, m)
% Calculates the results of TBT for a prolate ellipsoid by an
% interface.
%  INPUTS
% p.ep     -   slenderness parameter
% p.rho(s)   -   cross sectional shape
% p.rrhop(s) -   derivative of cross sectional shape times (rho(s)rho'(s))
% p.kappa(s) -   curvature of centerline
% p.r1(s)    -   x component of centerline
% p.r2(s)    -   y component of centerline
% p.r3(s)    -   z component of centerline
% p.t(s)     -   tangent to centerline
% p.erho1(s,phi) -   x component of local radial vector
% p.erho1(s,phi) -   y component of local radial vector
% p.erho1(s,phi) -   z component of local radial vector
% p.d        - distance below the interface
% p.lambda   - ratio of the viscosities of the two regions
% p.itau(s)  -   integrated torsion along the body
% p.p.tol - absolute p.tolerance for integrals
% m.cells - cells on which f is constant
% m.lookup - lookup for correct indices when indexing into m.cells
% m.numCells - number of cells, ie m.cells(:,1:m.numCells) is nontrivial.
% m.neighbours - cell of array of lists of neighbours for each cell.
%
% OUTPUTS
%  initMat        -   Matrix for initial solve.
%  iterMat        -   Matrix for subsequent iterations.

% Precompute coefficents needed for constructing the matrices.
lambdaB = (1 - p.lambda) / (1 + p.lambda);
lambda2 = 2 * p.lambda / (1 + p.lambda);
Amat = [1,0,0;0,1,0;0,0,-1];
Bmat = [lambdaB ,0,0;0,lambdaB,0;0,0,-1];

% Extract the bounds of the domain from the mesh.
sBounds = m.XBounds;
phiBounds = m.YBounds;

% Define the surface of the body.
X = @(s,phi) p.r1(s) + p.ep * p.rho(s) .* p.erho1(s,phi);
Y = @(s,phi) p.r2(s) + p.ep * p.rho(s) .* p.erho2(s,phi);
Z = @(s,phi) p.r3(s) + p.ep * p.rho(s) .* p.erho3(s,phi) - p.d;

% Construct the effective ellipsoids on each cell.
u = zeros(m.numCells,1);
c = zeros(m.numCells,1);
a = zeros(m.numCells,1);
alpha = zeros(m.numCells,1);
tv = zeros(m.numCells,3);
for cellInd = 1 : m.numCells
    s = m.cells(m.lookup.XMid, cellInd);
    phi = m.cells(m.lookup.YMid, cellInd);
    tv(cellInd,:) = p.t(s);
    rt = p.rho(s)^2 + sqrt(p.rho(s)^4 + 4*p.rrhop(s)^2);
    c(cellInd) = sqrt(rt/2);
    u(cellInd) =-2 * p.rrhop(s)/rt;
    a(cellInd) = 1 - p.ep * p.rho(s) * p.kappa(s) * cos(phi - p.itau(s));
    alpha(cellInd) = p.ep * c(cellInd) / a(cellInd);
end

% Compute zeta_1' and zeta_2' on each cell.
z1p = Z1p(a, alpha, u);
z2p = Z2p(a, alpha, u);

% Project the current 2D mesh onto a mesh for s. We'll use this to evaluate
% the phi-independent part of the matrix B.
mesh1D = projectMesh(m);

% Compute the s-dimensions of each cell in both the 1D and 2D mesh.
cellWidthS1D = mesh1D.cells(mesh1D.lookup.XWidth,:);
cellWidthS2D = m.cells(m.lookup.XWidth,:);
cellWidthPhi = m.cells(m.lookup.YWidth,:);

%% Compute the averages <1/zeta_1'> and <1/zeta_2'> over phi at each s in the
% 1D mesh.
z1pib = zeros(mesh1D.numCells,1);
z2pib = zeros(mesh1D.numCells,1);

% Form an interpolant over all the 2D cells, which will allow us to integrate
% even with non-uniform meshes.
int = scatteredInterpolant(reshape(m.cells(m.lookup.XMid,1:m.numCells),[],1), reshape(m.cells(m.lookup.YMid,1:m.numCells),[],1), 1./z1p, 'linear', 'nearest');
for sInd = 1 : mesh1D.numCells
    s = mesh1D.cells(mesh1D.lookup.XMid, sInd);
    z1pib(sInd) = integral(@(phi) int(s*ones(size(phi)),phi), phiBounds(1), phiBounds(2));
end
int.Values = 1./z2p;
for sInd = 1 : mesh1D.numCells
    s = mesh1D.cells(mesh1D.lookup.XMid, sInd);
    z2pib(sInd) = integral(@(phi) int(s*ones(size(phi)),phi), phiBounds(1), phiBounds(2));
end

%% Construct the matrix B. It is only a function of s.
BLocal = zeros(3*mesh1D.numCells); % This will be block diagonal.
BNonLocal = zeros(3*mesh1D.numCells); % This will be dense.

% Build the local and non-local parts of B, storing separately for later use.
textprogressbar('Building B:')
for sInd = 1 : mesh1D.numCells
    s = mesh1D.cells(mesh1D.lookup.XMid, sInd);
    % There is a local ellipsoidal contribution.
    % Form the outer product of the local centreline tangent with itself.
    T = p.t(s) .* p.t(s)';
    BLocal(3*(sInd-1) + (1:3), 3*(sInd-1) + (1:3)) = T / z1pib(sInd) + (eye(3) - T) / z2pib(sInd);

    % Compute the non-local contributions to all the 1D cells.
    for otherSInd = 1 : mesh1D.numCells
        % Perform the integrals of the various contributions of all the
        % kernels over a subinterval of the s domain. As the domain of
        % integration is shared, and it's 1D, we can compute all of these
        % integrals in one call to integral. Output is returned as a struct
        % with expected field names. If slow, one could use an ODE solver to
        % do this more efficiently.
        integrand = @(sOther) B_integrand(sOther, s, p.r1, p.r2, p.r3, p.rho, p.ep, p.d);
        o = unpack_B_integrals(integral(integrand, mesh1D.cells(mesh1D.lookup.XBoundLower, otherSInd), mesh1D.cells(mesh1D.lookup.XBoundUpper, otherSInd), 'ArrayValued', true, 'AbsTol', p.tol));

        BS = o.Bi * eye(3) +[o.B11 o.B12 o.B13; o.B12 o.B22 o.B23; o.B13 o.B23 o.B33];
        BSm =(o.Biw * eye(3) +[o.B11w o.B12w o.B13w; o.B12w o.B22w o.B23w; o.B13w o.B23w o.B33w]) * Bmat;
        BSDm = -lambda2 * (o.Biwsd * eye(3) -3 * [o.B11wsd o.B12wsd o.B13wsd; o.B12wsd o.B22wsd o.B23wsd; o.B13wsd o.B23wsd o.B33wsd]) * Amat;
        BSSm = -lambda2 * [0,0,o.B13wss;0,0,o.B23wss; -o.B13wss,-o.B23wss,0] * Amat;
        BNonLocal(((1:3)+(sInd-1) * 3),((1:3)+(otherSInd-1) * 3)) = BS + BSm + BSDm + BSSm;
    end
    textprogressbar(sInd / mesh1D.numCells * 100)
end
textprogressbar('.')
B = BLocal + BNonLocal;

%% Precompute local residuals that will be needed for builing the full matrix.
%% For each cell in the 2D mesh, we can compute its contribution to the cells
%% in the 1D mesh using B. Eventually, we'll add up contributions from the 1D
%% mesh.
MT0 = zeros(3,3,mesh1D.numCells,m.numCells);
residual = zeros(3*mesh1D.numCells, 3);
for cellInd = 1 : m.numCells

    % Compute quantities that are fixed on the cell.
    s = m.cells(m.lookup.XMid, cellInd);
    T = p.t(s) .* p.t(s)';
    Mi = cellWidthPhi(cellInd) * (T / z1p(cellInd) + (eye(3) - T) / z2p(cellInd));

    % B and BLocal are functions only of s, on the 1D mesh. We'll need to sum
    % up the contributions of any subcells of the 1D mesh that make up the
    % cell.
    residual(:) = 0;
    for subcellInd = mesh1D.projMap{cellInd}
        subcellEntryRange = (1:3) + 3*(subcellInd - 1);
        A0 = BLocal(:, subcellEntryRange) * Mi;
        residual = residual + (A0 - BLocal * (B \ A0));
    end

    % Save the precomputed residuals in a convenient format. Note that we've
    % computed the value for all the 1D cells at once (third index).
    MT0(:,1,:,cellInd) = reshape(residual(:,1),3,mesh1D.numCells);
    MT0(:,2,:,cellInd) = reshape(residual(:,2),3,mesh1D.numCells);
    MT0(:,3,:,cellInd) = reshape(residual(:,3),3,mesh1D.numCells);
end

% It is worth noting that the entries of MT0 do not scale with 1D mesh, ie
% taking a mesh with 2x as many elements does not reduce the size of the
% elements by approximately a factor of 2; instead, the element size is
% unchanged. Thus, we must take a weighted sum of their contributions below.

%% Build the matrix corresponding to the initialisation step, initMat. This
%% will use the precomputed values above, summing them over 1D subcells.
initMat = zeros(3*m.numCells);
textprogressbar('Building initMat:')
for cellInd = 1 : m.numCells
    cellEntryRange = (1:3) + 3*(cellInd - 1);
    s = m.cells(m.lookup.XMid, cellInd);
    T = p.t(s) .* p.t(s)';
    localContrib = (T / z1p(cellInd) + ((eye(3) - T) / z2p(cellInd)));

    % The cell may be made up of multiple 1D cells. We loop through the
    % relevant 1D cells and sum the contributions.
    for subcellInd = mesh1D.projMap{cellInd}
        % We add up the effects of the other cells on this subcell.
        subcellWeighting = cellWidthS1D(subcellInd) / cellWidthS2D(cellInd);
        for otherCellInd = 1 : m.numCells
            otherCellEntryRange = (1:3) + 3*(otherCellInd-1);
            initMat(cellEntryRange, otherCellEntryRange) = initMat(cellEntryRange, otherCellEntryRange) - localContrib * MT0(:,:,subcellInd,otherCellInd) * subcellWeighting;
        end
    end

    % Also insert the local contribution.
    initMat(cellEntryRange, cellEntryRange) = initMat(cellEntryRange, cellEntryRange) + localContrib;

    textprogressbar(cellInd / m.numCells * 100)
end
textprogressbar('.')

%% Compute the matrix iterMat needed for the iterative steps.
iterMat = zeros(3 * m.numCells);

% For each cell, we'll have to integrate the contributions of all the other
% cells. In order to parallelise this computation, we'll use more cumbersome
% linear indexing of cells and store the results.
% There are a total of m.numCells^2 ordered pairs.
numEntries = m.numCells^2;
% We'll store the results of the 23 necessary integrals per pair for later
% use.
computedIntegrals = zeros(23, numEntries);

% Parallel progress bar.
textprogressbar('Building iterMat:')
DQ = parallel.pool.DataQueue;
afterEach(DQ, @nUpdateProgressbar);
prog = 1;
parfor linIndex = 1 : numEntries
    % Get the indices of the cells from the linear index. One cell will define
    % the integral kernel (refCellInd), whilst the other defines the range of
    % integration (intCellInd).
    [refCellInd, intCellInd] = getCellIndices(linIndex, m.numCells);

    % To form the kernels efficiently, we precompute a range of quantities to
    % pass to the kernel function. These are all in terms of quantities on
    % refCellInd.

    % s and phi at the midpoint of the reference cell.
    sRef = m.cells(m.lookup.XMid, refCellInd);
    phiRef = m.cells(m.lookup.YMid, refCellInd);
    
    % Components of the centreline on the reference cell.
    r1 = p.r1(sRef);
    r2 = p.r2(sRef);
    r3 = p.r3(sRef);

    % Effective s point of the ellipsoid on the reference cell.
    uRef = u(refCellInd);
    r1e = @(sDummy) a(refCellInd) .* (sDummy-uRef) .* tv(refCellInd,1) + r1;
    r2e = @(sDummy) a(refCellInd) .* (sDummy-uRef) .* tv(refCellInd,2) + r2;
    r3e = @(sDummy) a(refCellInd) .* (sDummy-uRef) .* tv(refCellInd,3) + r3;

    % Radius of the ellipsoid at s = sDummy.
    rhoEll = @(sDummy) c(refCellInd) .* sqrt(1-sDummy.^2);
    % Surface points on the ellipsoid.
    XEll = @(sDummy,phiDummy) r1e(sDummy) + p.ep * rhoEll(sDummy) .* p.erho1(sRef,phiDummy);
    YEll = @(sDummy,phiDummy) r2e(sDummy) + p.ep * rhoEll(sDummy) .* p.erho2(sRef,phiDummy);
    ZEll = @(sDummy,phiDummy) r3e(sDummy) + p.ep * rhoEll(sDummy) .* p.erho3(sRef,phiDummy);
    sEll = @(sDummy) uRef - sRef + sDummy;

    % Define the domain of integration.
    lowSCell = m.cells(m.lookup.XBoundLower, intCellInd);
    uppSCell = m.cells(m.lookup.XBoundUpper, intCellInd);
    lowPhiCell = m.cells(m.lookup.YBoundLower, intCellInd);
    uppPhiCell = m.cells(m.lookup.YBoundUpper, intCellInd);
    
    % Free-space contributions.
    % Define intermediate, per-thread storage.
    toStore = zeros(23, 1);
    C = zeros(7,1);
    if refCellInd == intCellInd 

        % If we are integrating over the reference cell, there is a singularity in
        % the domain of integration. Hence, split up the domain.

        % Compute additional upper and lower bounds for the domain of
        % integration, noting that some integrals are only over the surface
        % of the ellipsoid.
        lowSEll = -1 + sRef - uRef;
        lowSThresh = max(lowSCell, lowSEll);

        uppSEll = 1 + sRef - uRef;
        uppSThresh = min(uppSCell, uppSEll);
    
        C = C + quad2dv(@(s,phi) c1Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),lowSThresh,sRef,lowPhiCell,uppPhiCell,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c2Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),sRef,uppSThresh,lowPhiCell,uppPhiCell,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c3Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),lowSEll,lowSThresh,phiBounds(1),phiBounds(2),'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c4Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),lowSCell,lowSThresh,lowPhiCell,uppPhiCell,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c5Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),uppSThresh,uppSEll,phiBounds(1),phiBounds(2),'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c6Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),uppSThresh,uppSCell,lowPhiCell,uppPhiCell,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c7Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),lowSThresh,uppSThresh,phiBounds(1),lowPhiCell,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c8Integrand(s,phi,sRef,phiRef,uRef,sEll,X,Y,Z,XEll,YEll,ZEll),lowSThresh,uppSThresh,uppPhiCell,phiBounds(2),'AbsTol',p.tol);

    else

        % The integration and reference cells are distinct, so the integrand is regular.
        integrand = @(x,y) cStokesletIntegrand(x,y,sRef,phiRef,X,Y,Z,XEll,YEll,ZEll);
        C = quad2dv(integrand,lowSCell,uppSCell,lowPhiCell,uppPhiCell,'AbsTol',p.tol);

    end
    
    % Contributions due to the boundary.
    integrand = @(x,y) CWallIntegrand(x,y,sRef,phiRef,X,Y,Z,XEll,YEll,ZEll);
    D = quad2dv(integrand,lowSCell,uppSCell,lowPhiCell,uppPhiCell,'AbsTol',p.tol); 

    % Store the computed integrals. This data structure is cumbersome, but
    % efficient. The indices are as follows, in terms of original variable
    % names:
    % 1  : ci
    % 2  : c11
    % 3  : c12
    % 4  : c13
    % 5  : c22
    % 6  : c23
    % 7  : c33
    % 8  : ciw
    % 9  : c11w
    % 10 : c12w
    % 11 : c13w
    % 12 : c22w
    % 13 : c23w
    % 14 : c33w
    % 15 : ciwsd
    % 16 : c11wsd
    % 17 : c12wsd
    % 18 : c13wsd
    % 19 : c22wsd
    % 20 : c23wsd
    % 21 : c33wsd
    % 22 : c13wss
    % 23 : c23wss

    toStore(1:7) = C(1:7);
    toStore(8:23) = D(1:16);
    computedIntegrals(:,linIndex) = toStore;

    send(DQ, linIndex);
end
textprogressbar('.')

% Use the computed values to construct the actual matrix. We'll do this via an
% intermediate, redefining C from above.
C = zeros(3 * m.numCells);
for linIndex = 1 : numEntries
    % Get the indices of the cells from the linear index. One cell will define
    % the integral kernel (refCellInd), whilst the other defines the range of
    % integration (intCellInd).
    [refCellInd, intCellInd] = getCellIndices(linIndex, m.numCells);
    cI = computedIntegrals(:,linIndex);

    % Construct some further intermediates.
    CS = cI(1) * eye(3) + symmMat(cI(2:7),3);
    CSm = (cI(8) * eye(3) + symmMat(cI(9:14),3)) * Bmat;
    CSDm = - lambda2 * (cI(15) * eye(3) - 3 * symmMat(cI(16:21),3)) * Amat;
    CSSm = - lambda2 * [0,0,cI(22); 0,0,cI(23); -cI(22), -cI(23),0] * Amat;

    % We need to sum up the contributions of various 1D subcells from
    % BNonlocal. We need to add up the subcells of both the reference cell
    % and the integration cell, weighting the contributions of the
    % integration subcells.
    BNonlocalSum = zeros(3);
    for refSubcellInd = mesh1D.projMap{refCellInd}
        refSubcellEntryRange = (1:3) + 3*(refSubcellInd - 1);

        for intSubcellInd = mesh1D.projMap{intCellInd}
            intSubcellEntryRange = (1:3) + 3*(intSubcellInd - 1);
            % Define the weighting factor of the integration subcell.
            subcellWeighting = cellWidthS1D(intSubcellInd) / cellWidthS2D(intCellInd);
            BNonlocalSum = BNonlocalSum + BNonLocal(refSubcellEntryRange, intSubcellEntryRange) * subcellWeighting;
        end
        
    end

    refCellEntryRange = (1:3) + 3*(refCellInd - 1);
    intCellEntryRange = (1:3) + 3*(intCellInd - 1);
    C(refCellEntryRange, intCellEntryRange) = CS - BNonlocalSum * cellWidthPhi(intCellInd) + CSm + CSDm + CSSm;

    % If the reference and integration cells are the same, add in a local contribution.
    if refCellInd == intCellInd
        sRef = m.cells(m.lookup.XMid, refCellInd);
        T = p.t(sRef) .* p.t(sRef)';
        localContrib = T * (z1p(refCellInd) - Z1(a(refCellInd),alpha(refCellInd))) + (eye(3) - T) * (z2p(refCellInd) - Z2(a(refCellInd),alpha(refCellInd))); 
        C(refCellEntryRange, intCellEntryRange) = C(refCellEntryRange, intCellEntryRange) - localContrib;
    end

end

% Form the iterMat from the above.
iterMat = initMat * C;

function nUpdateProgressbar(~)
        textprogressbar( prog / numEntries * 100);
        prog = prog + 1;
end

end