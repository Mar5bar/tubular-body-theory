addpath(genpath('helpers'))
%% Setup

% Maximum number of mesh refinements.
maxRefs = 5;

% Threshold for refinement. Lower leads to more refinement.
refThresh = 0.25;

% Slenderness of the body.
ep = 1;

% Cross-sectional radius; function of arclength.
rho = @(s) (1-s.^2).^(0.5);

% rho * d/ds (rho); function of arclength.
rrhop = @(s) -s;

% Centreline curvature; function of arclength.
kappa = @(s) 0*s;

% Cartesian components of the centreline; function of arclength.
r1 = @(s) s;
r2 = @(s) 0*s;
r3 = @(s) 0*s;

% Tangent vector to the centreline; function of arclength.
t = @(s) [1;0;0] + 0*s;

% Cartesian components of the local radial vector; function of arclength and angle.
erho1 = @(s,phi) 0*s;
erho2 = @(s,phi) cos(phi);
erho3 = @(s,phi) sin(phi);

% Integral of the torsion along the centreline.
itau = @(s) 0*s;

% z offset from the plane boundary.
d = 1.1;

% Ratio of viscosity of the two fluid regions.
lambda = 1;

% Number of iterations when computing f (>=0).
numIter = 1000;

% Absolute tolerance for integrals.
tol = 1e-6;

% Assign parameters to a struct for easy passing.
params.ep = ep;
params.rho = rho;
params.rrhop = rrhop;
params.kappa = kappa;
params.r1 = r1;
params.r2 = r2;
params.r3 = r3;
params.t = t;
params.erho1 = erho1;
params.erho2 = erho2;
params.erho3 = erho3;
params.itau = itau;
params.d = d;
params.lambda = lambda;
params.tol = tol;
params.numIter = numIter;
params.refThresh = refThresh;

%% Generate an initial mesh of cells.
% Bounds on S and Phi.
sBounds = [-1,1];
phiBounds = [0,2*pi];

% Initial number of cells in the S and Phi coordinates..
initSNum = 5;
initPhiNum = 4;

% Number of cells to allocate initially.
cellsToAllocate = 1e3;
assert(initSNum * initPhiNum < cellsToAllocate, 'Initial discretisation exceeds cellsToAllocate.')

% Generate the mesh.
mesh2D = initMesh(sBounds, phiBounds, initSNum, initPhiNum, cellsToAllocate);

%% Compute the resistance matrices.
% We'll iteratively refine the coarse mesh to generate smooth approximations to f.
numRefs = -1;
refMask = [];
meshes = {};
meshes{1} = mesh2D;
Rs = {};
fTotals = {};
while (numRefs < maxRefs && any(refMask)) || numRefs == -1

    % If this isn't the first pass, refine the mesh.
    if numRefs >= 0
        mesh2D = refineMesh(mesh2D, refInds);
        meshes{numRefs+2} = mesh2D;
    end

    % Compute the matrices needed for finding f.
    [initMat, iterMat] = TBT_interface(params,mesh2D);
    [R,f,fTotal] = Rmat(params,mesh2D,initMat,iterMat);
    Rs{numRefs+2} = R(:,:,end);
    fTotals{numRefs+2} = fTotal{end};

    % Test, via fsTotal, if we need to refine the mesh.
    [refMask, refInds] = refinementNeeded(fTotal{end}, mesh2D.neighbours, params.refThresh);

    numRefs = numRefs + 1;

end