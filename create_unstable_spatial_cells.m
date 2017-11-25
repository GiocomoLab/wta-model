function [maps] = create_unstable_spatial_cells(smallest_unstable_cell,peakFR_all,numUnstableSpatialcells)

%% create the grid patterns
mapSize = 100; boxSize = 50;
fieldSize = ceil(2000*rand(numUnstableSpatialcells,1))+smallest_unstable_cell;
peakFR = peakFR_all(ceil(length(peakFR_all)*rand(numUnstableSpatialcells,1)));
[X1,X2] = meshgrid(linspace(0,boxSize,mapSize),linspace(0,boxSize,mapSize));
muVec = [ceil(2*boxSize*rand(numUnstableSpatialcells,1))-boxSize/2 ceil(2*boxSize*rand(numUnstableSpatialcells,1))-boxSize/2];

% initialize matrices

maps = zeros(numUnstableSpatialcells,mapSize^2);
% create Gaussian bumps "grid cells"
for k = 1:numUnstableSpatialcells
    mu = muVec(k,:);
    Sigma = fieldSize(k)/(-2*log(0.2)*pi)*eye(2);
    field = mvnpdf([X1(:) X2(:)],mu,Sigma);
    pattern = peakFR(k).*field./max(max(field));
    maps(k,:) = pattern;
end