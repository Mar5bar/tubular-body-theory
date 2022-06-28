addpath(genpath('helpers'))
f = @(x,y) 1./((x.^2 + y.^2 + 0.001));

% Maximum allowed number of cells.
maxCells = 10000;

% Threshold for refinement.
thresh = 0.25;

% Initial discretisation parameters.
initXNum = 10;
initYNum = 10;
assert(initXNum * initYNum < maxCells, 'Initial discretisation exceeds maxCells.')

bounds = struct();
bounds.X = [-1,1];
bounds.Y = [-1,1];

[cells, lookup, numCells, neighbours] = initMesh(bounds, initXNum, initYNum, maxCells);

fs = zeros(numCells,1);
for cellInd = 1 : numCells
    fs(cellInd) = f(cells(lookup.XMid,cellInd), cells(lookup.YMid,cellInd));
end

[refMask, refInds] = refinementNeeded(fs, neighbours, thresh);

figure
nexttile()
hold on
scatter3(cells(lookup.XMid,1:numCells), cells(lookup.YMid,1:numCells), fs(1:numCells), [], fs(1:numCells), 'filled')
scatter3(cells(lookup.XMid,refInds), cells(lookup.YMid,refInds), fs(refInds), [], fs(refInds), 'MarkerEdgeColor','black','LineWidth',3)
colormap(viridis)
view(0,90)
drawCells(cells,lookup)

iterCount = 0;
while any(refMask) & numCells <= maxCells

    [cells, numCells, neighbours] = refineCells(cells, lookup, numCells, refInds);

    fs = zeros(numCells,1);
    for cellInd = 1 : numCells
        fs(cellInd) = f(cells(lookup.XMid,cellInd), cells(lookup.YMid,cellInd));
    end

    [refMask, refInds] = refinementNeeded(fs, neighbours, thresh);

    iterCount = iterCount + 1;

    nexttile()
    hold on
    scatter3(cells(lookup.XMid,1:numCells), cells(lookup.YMid,1:numCells), fs(1:numCells), [], fs(1:numCells), 'filled')
    scatter3(cells(lookup.XMid,refInds), cells(lookup.YMid,refInds), fs(refInds), [], fs(refInds), 'MarkerEdgeColor','black','LineWidth',3)
    colormap(viridis)
    view(0,90)
    drawCells(cells,lookup)

end

if numCells > maxCells
    disp('Maximum number of cells exceeded')
end
if ~any(refMask)
    disp(['Variation within tolerance after ',num2str(iterCount),' iterations'])
end