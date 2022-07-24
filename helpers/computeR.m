function R = computeR(f, m, pos)
%% computeR will form the resistance matrix corresponding to the vector of
%% forces f on a mesh m, with cell centres in pos.
% First, weight the forces (per unit area) by the area of each cell.
weightedF = f .* kron(m.cells(m.lookup.XWidth,1:m.numCells)' .* m.cells(m.lookup.YWidth,1:m.numCells)', [1;1;1]);
R(1:3,:,1) = sum(reshape(weightedF, 3, [], 6), 2); % Total force.
R(4:6,:,1) = sum(cross(repmat(pos,1,1,6),reshape(weightedF,3,[],6)),2); % Total torque.