addpath(genpath('helpers'))
%% Setup

% Maximum number of mesh refinements.
maxRefs = 5;

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

% Number of iterations (>=0).
N = 1000;

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

%% Generate an initial mesh of cells.
% Bounds on S and Phi.
sBounds = [-1,1];
phiBounds = [0,2*pi];

% Initial number of cells in the S and Phi coordinates..
initSNum = 5;
initPhiNum = 5;

% Number of cells to allocate initially.
cellsToAllocate = 1e3;
assert(initSNum * initPhiNum < cellsToAllocate, 'Initial discretisation exceeds cellsToAllocate.')

% Generate the mesh.
mesh2D = initMesh(sBounds, phiBounds, initSNum, initPhiNum, cellsToAllocate);

% Perform the initial computation on the coarse mesh.
% Compute the matrices needed for finding f.
[initMat,iterMat] = TBT_interface(params,mesh2D);
[R,fs,fsTotal,S,Phi] = Rmat(initMat,iterMat,N,numS,numPhi,epps,rho,r1,r2,r3,erho1,erho2,erho3);

% We'll iteratively refine the coarse mesh to generate smooth approximations to f.
% The first iteration will compute the solution and will not refine.
numRefs = 0;
while (numRefs < maxRefs & any(refMask)) | numRefs == 0

    % Compute the matrices needed for finding f.
    [SO0, SNO] = TBT_interface(params,mesh2D);
    [R,fs,fsTotal,S,Phi] = Rmat(initMat,iterMat,N,numS,numPhi,epps,rho,r1,r2,r3,erho1,erho2,erho3);

    numRefs = numRefs + 1;
end




%% Call TBT_interface
tic
[SO0,SNO] = TBT_interface(params, mesh2D);
toc
% tic
% [SO0OLD,SNOOLD] = TBT_interfacev3OLD(epps,rho,rrhop,kappa,r1,r2,r3,t,erho1,erho2,erho3,itau,d, lambda, N, numS);
% toc

disp('Iterating')
tic
[R,fs,fsTotal,S,Phi] = Rmat(SO0,SNO,N,numS,numPhi,epps,rho,r1,r2,r3,erho1,erho2,erho3);
toc
S = S'; Phi = Phi';

% Plot the final approximation to f in the 6 test cases needed to compute R.
plotF(fsTotal{end},S,Phi)