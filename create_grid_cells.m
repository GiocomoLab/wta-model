function [allMaps_grid,mapSize] = create_grid_cells(smallest_grid,peakFR_all,numGridCells)

%% create the grid patterns

% set the parameters for the grid cells
mapSize = 100; boxSize = 50;
peakFR = peakFR_all(ceil(length(peakFR_all)*rand(numGridCells,1)));
allMaps_grid = nan(numGridCells,mapSize^2);
orientation = 7*ones(numGridCells,1);
gridSpacing = [smallest_grid smallest_grid*sqrt(2) smallest_grid*2 smallest_grid*2*sqrt(2) smallest_grid*4]';
spacingVec = gridSpacing(ceil(length(gridSpacing)*rand(numGridCells,1)));
fieldSizeVec = pi*(0.22*spacingVec+3.95).^2;
phaseVec = rand(numGridCells,1);

% create grid pattern for each cell
for k = 1:numGridCells
    gridMap  = create_grid_pattern(orientation(k),spacingVec(k),phaseVec(k),boxSize,peakFR(k),fieldSizeVec(k));
    allMaps_grid(k,:) = gridMap(:)';
end

