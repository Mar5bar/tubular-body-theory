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
cellWidthS1D = mesh1D.cells(mesh1D.lookup.XBoundUpper,:) - mesh1D.cells(mesh1D.lookup.XBoundLower,:);
cellWidthS2D = m.cells(m.lookup.YBoundUpper,:) - m.cells(m.lookup.YBoundLower,:);

%% Compute the averages <1/zeta_1'> and <1/zeta_2'> over phi at each s in the
% 1D mesh.
z1pib = zeros(mesh1D.numCells,1);
z2pib = zeros(mesh1D.numCells,1);

% Form an interpolant over all the 2D cells, which will allow us to integrate
% even with non-uniform meshes.
int = scatteredInterpolant(reshape(m.cells(m.lookup.XMid,1:m.numCells),[],1), reshape(m.cells(m.lookup.YMid,1:m.numCells),[],1), 1./z1p);
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
    Mi = T / z1p(cellInd) + (eye(3) - T) / z2p(cellInd);

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

end


% -------------------
numS = 5;
numPhi = 5;
dphi = 2*pi/5;
Mi = zeros(3,3,numPhi,numS,numPhi,numS);
for j = 1:numS
    for i = 1:numPhi
        Mi(:,:,i,j,i,j) = (p.t(s(j)) .* p.t(s(j))')/z1p(i+5*(j-1)) +((eye(3)-p.t(s(j)) .* p.t(s(j))'))/z2p(i+5*(j-1));
    end
end

MT0 = zeros(3,3,numS,numPhi,numS);
for k = 1:numS
    for l2 = 1:numPhi
        Q1 = reshape(dphi * sum(Mi(:,1,:,:,l2,k),3),3 * numS,1);
        Q2 = reshape(dphi * sum(Mi(:,2,:,:,l2,k),3),3 * numS,1);
        Q3 = reshape(dphi * sum(Mi(:,3,:,:,l2,k),3),3 * numS,1);
        
        A0 = BLocal * [Q1, Q2, Q3];

        fb = B \ A0;
        
        res = A0 - BLocal * fb;
        MT0(:,1,:,l2,k) = reshape(res(:,1),3,numS);
        MT0(:,2,:,l2,k) = reshape(res(:,2),3,numS);
        MT0(:,3,:,l2,k) = reshape(res(:,3),3,numS);

    end
end

% Construct initMat from the above.
initMat = zeros(3*m.numCells);
textprogressbar('Building initMat:')
for k = 1:numS
    for l2 = 1:numPhi % Loop through all cells.
        for j = 1 : numS
            T = t(s(j)) .* t(s(j))';
            for i = 1:numPhi
                MiContrib = (T/z1p(i,j) +((eye(3)-T)/z2p(i,j)));
                initMat((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = MiContrib*(i==l2 && j==k) - MiContrib * MT0(:,:,j,l2,k);
            end
        end
    end
    textprogressbar(k / numS * 100)
end
textprogressbar('.')


%% Next order creation
    
%Create the full matrix for integral equations
ci  = zeros(numPhi * numS * numPhi * numS,1);
c11 = zeros(numPhi * numS * numPhi * numS,1);
c12 = zeros(numPhi * numS * numPhi * numS,1);
c13 = zeros(numPhi * numS * numPhi * numS,1);
c22 = zeros(numPhi * numS * numPhi * numS,1);
c23 = zeros(numPhi * numS * numPhi * numS,1);
c33 = zeros(numPhi * numS * numPhi * numS,1);

ciw  = zeros(numPhi * numS * numPhi * numS,1);
c11w = zeros(numPhi * numS * numPhi * numS,1);
c12w = zeros(numPhi * numS * numPhi * numS,1);
c13w = zeros(numPhi * numS * numPhi * numS,1);
c22w = zeros(numPhi * numS * numPhi * numS,1);
c23w = zeros(numPhi * numS * numPhi * numS,1);
c33w = zeros(numPhi * numS * numPhi * numS,1);

ciwsd  = zeros(numPhi * numS * numPhi * numS,1);
c11wsd = zeros(numPhi * numS * numPhi * numS,1);
c12wsd = zeros(numPhi * numS * numPhi * numS,1);
c13wsd = zeros(numPhi * numS * numPhi * numS,1);
c22wsd = zeros(numPhi * numS * numPhi * numS,1);
c23wsd = zeros(numPhi * numS * numPhi * numS,1);
c33wsd = zeros(numPhi * numS * numPhi * numS,1);

c13wss = zeros(numPhi * numS * numPhi * numS,1);
c23wss = zeros(numPhi * numS * numPhi * numS,1);

numEntries = numPhi * numS * numPhi * numS;
s = linspace(-1+(ds/2),1-(ds/2),numS);
phi = linspace(-pi+dphi/2,pi-dphi/2,numPhi);

DQ = parallel.pool.DataQueue;
textprogressbar('Building iterMat:')
afterEach(DQ, @nUpdateProgressbar);
prog = 1;
parfor l = 1 : numEntries 
    [k,l2,j,i,sameCell] = getIndex(l,numS,numPhi);
    
    pt = phi(i);
    st = s(j);
    
    a2 = a(i,j);
    t12 = tv(i,j,1);
    t22 = tv(i,j,2);
    t32 = tv(i,j,3);
    r12 = r1(st);
    r22 = r2(st);
    r32 = r3(st);
    c2 = c(i,j);
    u2 = u(i,j); %effective s point of ellipsoid
    r1e = @(s2) a2 .* (s2-u2) .* t12 +r12;
    r2e = @(s2) a2 .* (s2-u2) .* t22 +r22;
    r3e = @(s2) a2 .* (s2-u2) .* t32 +r32;
    rhoe = @(s2) c2 .* sqrt(1-s2.^2);
    Xe = @(s2,phi2) r1e(s2) + ep * rhoe(s2) .* erho1(st,phi2);
    Ye = @(s2,phi2) r2e(s2) + ep * rhoe(s2) .* erho2(st,phi2);
    Ze = @(s2,phi2) r3e(s2) + ep * rhoe(s2) .* erho3(st,phi2);
    s3 = @(s2) u2-st+s2;
    
    % Normal contributions
    C = zeros(7,1);
    if sameCell % Cell with the singularity within it
        
        low = max(s(k)-ds/2,-1+st-u2);
        upp = min(s(k)+ds/2,1+st-u2);
    
        C = C + quad2dv(@(s,phi) c1Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),low,st,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c2Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),st,upp,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c3Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),-1+st-u2,low,-pi,pi,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c4Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),s(k)-ds/2,low,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c5Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),upp,1+st-u2,-pi,pi,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c6Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),upp,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c7Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),low,upp,-pi,phi(l2)-dphi/2,'AbsTol',p.tol);
        C = C + quad2dv(@(s,phi) c8Integrand(s,phi,st,pt,u2,s3,X,Y,Z,Xe,Ye,Ze),low,upp,phi(l2)+dphi/2,pi,'AbsTol',p.tol);

    else % If there is no singularity within the cell we can evaluate the terms directly

        integrand = @(x,y) cStokesletIntegrand(x,y,st,pt,X,Y,Z,Xe,Ye,Ze);
        C = quad2dv(integrand,s(k)-ds/2,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol);
        
    end

    ci(l)  = C(1);
    c11(l) = C(2);
    c12(l) = C(3);
    c13(l) = C(4);
    c22(l) = C(5);
    c23(l) = C(6);
    c33(l) = C(7);
    
    % INTERFACE CONTRIBUTIONS
    integrand = @(x,y) CWallIntegrand(x,y,st,pt,X,Y,Z,Xe,Ye,Ze);
    D = quad2dv(integrand,s(k)-ds/2,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',p.tol); 
    ciw(l)    = D(1);
    c11w(l)   = D(2);
    c12w(l)   = D(3);
    c13w(l)   = D(4);
    c22w(l)   = D(5);
    c23w(l)   = D(6);
    c33w(l)   = D(7);
    ciwsd(l)  = D(8);
    c11wsd(l) = D(9);
    c12wsd(l) = D(10);
    c13wsd(l) = D(11);
    c22wsd(l) = D(12);
    c23wsd(l) = D(13);
    c33wsd(l) = D(14);
    c13wss(l) = D(15);
    c23wss(l) = D(16);
    
    send(DQ, l);
end
textprogressbar('.')

%% CONSTRUCTING MATRIX SYSTEM.
C = zeros(3 * numS * numPhi,3 * numS * numPhi); %Construct the matrix C

for l = 1 : numEntries
    [k,l2,j,i,sameCell] = getIndex(l, numS, numPhi);
    
    CS =ci(l) * eye(3) +[c11(l) c12(l) c13(l); c12(l) c22(l) c23(l); c13(l) c23(l) c33(l)];
    CSm =(ciw(l) * eye(3) +[c11w(l) c12w(l) c13w(l); c12w(l) c22w(l) c23w(l); c13w(l) c23w(l) c33w(l)]) * Bmat;
    CSDm = -lambda2 * (ciwsd(l) * eye(3) -3 * [c11wsd(l) c12wsd(l) c13wsd(l); c12wsd(l) c22wsd(l) c23wsd(l); c13wsd(l) c23wsd(l) c33wsd(l)]) * Amat;
    CSSm = -lambda2 * [0,0,c13wss(l);0,0,c23wss(l); -c13wss(l),-c23wss(l),0] * Amat;
% BNonlocal will need its indexing adjusted - it current indexes into s values, but these need mapping via sUniqueMid.
    C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = CS -BNonlocal(((1:3)+(j-1) * 3),((1:3)+(k-1) * 3)) * dphi+CSm+CSDm+CSSm;
    if sameCell
        C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) -( (t(s(j)) .* t(s(j))') * (z1p(i,j)-Z1(a(i,j),alpha(i,j))) +((eye(3)-t(s(j)) .* t(s(j))')) * (z2p(i,j)-Z2(a(i,j),alpha(i,j))));
    end
end

% Form the iterMat from the above.
iterMat = initMat * C;

function nUpdateProgressbar(~)
        textprogressbar( prog / numEntries * 100);
        prog = prog + 1;
end

end