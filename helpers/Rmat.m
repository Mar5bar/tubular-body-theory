function [R, fs, fsTotal] = Rmat(p, m, initMat, iterMat)
%% Compute the resistance matrix using initMat and iterMat. Parameters are
%% passed in the structure p, and the mesh associated with initMat and
%% iterMat is stored in the 2D mesh m.

% We'll form and store the 6x6 resistance matrices computed at each
% iteration.
R = zeros(6, 6, p.numIter + 1);

% Define the surface of the body.
X = @(s,phi) p.r1(s) + p.ep * p.rho(s) .* p.erho1(s,phi);
Y = @(s,phi) p.r2(s) + p.ep * p.rho(s) .* p.erho2(s,phi);
Z = @(s,phi) p.r3(s) + p.ep * p.rho(s) .* p.erho3(s,phi);


% The evaluation points (s,phi) are as stored in the mesh.
s = m.cells(m.lookup.XMid,1:m.numCells);
phi = m.cells(m.lookup.YMid,1:m.numCells);
% Compute and store the positions to allow for fast, matrix-based operations below.
pos = [X(s,phi); Y(s,phi); Z(s,phi)];

% Define the surface velocity in the six cases, stored as columns of a matrix.
V = zeros(3 * m.numCells, 6);
V(:,1) = repmat([1;0;0], numel(s), 1);
V(:,2) = repmat([0;1;0], numel(s), 1);
V(:,3) = repmat([0;0;1], numel(s), 1);
V(:,4) = reshape(cross(repmat([1;0;0],1,numel(s)), pos), [], 1);
V(:,5) = reshape(cross(repmat([0;1;0],1,numel(s)), pos), [], 1);
V(:,6) = reshape(cross(repmat([0;0;1],1,numel(s)), pos), [], 1);
V = 8*pi * V;

% Compute the initial solution for f.
fs = cell(p.numIter + 1, 1);
fsTotal = cell(p.numIter + 1, 1);

fs{1} = initMat * V;
fsTotal{1} = fs{1};

% Compute a resistance matrix R(:,:,1) using the first approximation.
% First, weight the forces (per unit area) by the area of each cell.
R(:,:,1) = computeR(fs{1}, m, pos);

% Iterate to construct R using f_n.
sgn = -1;
for n = 1 : p.numIter
    fs{n+1} = iterMat * fs{n};
    fsTotal{n+1} = fsTotal{n} + sgn * fs{n+1};
    R(:,:,n+1) = computeR(fsTotal{n+1}, m, pos);
    sgn = -sgn;
end

end

