addpath(genpath('helpers'))
f = @(x,y) 1./((x.^2 + y.^2 + 0.001));

% Maximum allowed number of m.cells.
maxCells = 10000;

% Threshold for refinement.
thresh = 0.25;

% Initial discretisation parameters.
initXNum = 10;
initYNum = 10;
assert(initXNum * initYNum < maxCells, 'Initial discretisation exceeds maxCells.')

XBounds = [-1,1];
YBounds = [-1,1];

m = initMesh(XBounds, YBounds, initXNum, initYNum, maxCells);

fs = zeros(m.numCells,1);
for cellInd = 1 : m.numCells
    fs(cellInd) = f(m.cells(m.lookup.XMid,cellInd), m.cells(m.lookup.YMid,cellInd));
end

[refMask, refInds] = refinementNeeded(fs, m.neighbours, thresh);

figure
nexttile()
hold on
scatter3(m.cells(m.lookup.XMid,1:m.numCells), m.cells(m.lookup.YMid,1:m.numCells), fs(1:m.numCells), [], fs(1:m.numCells), 'filled')
scatter3(m.cells(m.lookup.XMid,refInds), m.cells(m.lookup.YMid,refInds), fs(refInds), [], fs(refInds), 'MarkerEdgeColor','black','LineWidth',3)
colormap(viridis)
view(0,90)
drawCells(m.cells,m.lookup)

iterCount = 0;
while any(refMask) & m.numCells <= maxCells

    [m.cells, m.numCells, m.neighbours] = refineCells(m.cells, m.lookup, m.numCells, refInds);

    fs = zeros(m.numCells,1);
    for cellInd = 1 : m.numCells
        fs(cellInd) = f(m.cells(m.lookup.XMid,cellInd), m.cells(m.lookup.YMid,cellInd));
    end

    [refMask, refInds] = refinementNeeded(fs, m.neighbours, thresh);

    iterCount = iterCount + 1;

    nexttile()
    hold on
    scatter3(m.cells(m.lookup.XMid,1:m.numCells), m.cells(m.lookup.YMid,1:m.numCells), fs(1:m.numCells), [], fs(1:m.numCells), 'filled')
    scatter3(m.cells(m.lookup.XMid,refInds), m.cells(m.lookup.YMid,refInds), fs(refInds), [], fs(refInds), 'MarkerEdgeColor','black','LineWidth',3)
    colormap(viridis)
    view(0,90)
    drawCells(m.cells,m.lookup)

end

if m.numCells > maxCells
    disp('Maximum number of m.cells exceeded')
end
if ~any(refMask)
    disp(['Variation within tolerance after ',num2str(iterCount),' iterations'])
end