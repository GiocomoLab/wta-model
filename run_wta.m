%% clear the workspace
clear all; close all; clc

%% load peak firing rate distribution
% load the peak firing rate distribution used for both populations
% This data were taken directly from grid cells recorded from Mallory et al,
% 2017. Cells with a firing rate > 10 Hz were ignored.
load peakFR_all;

%% set parameters for model

% smallest_grid = size of the smallest grid scale module
% smallest_unstable_spatial = size of smallest unstable spatial cells (cm^2)
% num_grid = # of grid cell projections that each place cell receives
% num_unstable_spatial = # of unstable spatial cell projections that each place cell receives
% numPlace = # of place cells in the model
% numGridCells = # of grid cells (total) in the model
% numUnstableSpatialcells = # of unstable spatial cells (total) in the model

smallest_grid = 25:5:60;
grid_scale_iter = numel(smallest_grid);
smallest_unstable_spatial = 1000;
num_grid = 1000;
num_unstable_spatial = 200;
numPlaceCells = 2000;
numGridCells = 3000;
numUnstableSpatialcells = 5000;
mapSize = 100;
E = 0.1;

%% initialize matrices
corr_bw_days = nan(numPlaceCells,grid_scale_iter);
size_place1 = nan(numPlaceCells,grid_scale_iter);
size_place2 = nan(numPlaceCells,grid_scale_iter);

%% run the WTA model!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE UNSTABLE SPATIAL INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the LEC maps
unstable_spatial_maps = create_unstable_spatial_cells(smallest_unstable_spatial,peakFR_all,numUnstableSpatialcells);

% compute the unstable spatial cell --> place cell connection probability
us_prob = num_unstable_spatial/numUnstableSpatialcells;

% compute the weight matrix between unstable spatial cells and place cells for day 1
W_usp1 = rand(numPlaceCells,numUnstableSpatialcells);
W_usp1(rand(numPlaceCells,numUnstableSpatialcells)>us_prob) = 0;

% and again for day 2
W_usp2 = rand(numPlaceCells,numUnstableSpatialcells);
W_usp2(rand(numPlaceCells,numUnstableSpatialcells)>us_prob) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE GRID WEIGHT MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the grid --> place connection probability and weight matrix
gp_prob = num_grid/numGridCells;
W_gp = rand(numPlaceCells,numGridCells);
W_gp(rand(numPlaceCells,numGridCells)>gp_prob) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN WTA MODEL FOR EACH SCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:grid_scale_iter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE GRID INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create the grid maps (used for both days)
    grid_maps = create_grid_cells(smallest_grid(n),peakFR_all,numGridCells);
    
    grid_input = W_gp*grid_maps;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIRST SESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    unstable_input1 = W_usp1*unstable_spatial_maps;
    
    day1_input = grid_input + unstable_input1;
    
    [day1_output] = wta_mechanism(day1_input,E);
    
    % compute the size of largest place field for each place cell
    for j = 1:numPlaceCells
        if any(day1_output(j,:) > 0)
            stats1 = regionprops(reshape(day1_output(j,:),mapSize,mapSize)>0,'Area');
            area1 = cat(1, stats1.Area);
            size_place1(j,n) = max(area1);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECOND SESSION %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    unstable_input2 = W_usp2*unstable_spatial_maps;
    
    day2_input = grid_input + W_usp2*unstable_spatial_maps;
    
    [day2_output] = wta_mechanism(day2_input,E);
    
    % compute the size of largest place field for each place cell
    for j = 1:numPlaceCells
        if any(day2_output(j,:) > 0)
            stats2 = regionprops(reshape(day2_output(j,:),mapSize,mapSize)>0,'Area');
            area2 = cat(1, stats2.Area);
            size_place2(j,n) = max(area2);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE THE TWO %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    corr_bw_days(:,n) = diag(corr(day1_output',day2_output'));

end

%% Basic analysis
% take away sessions with fields that are too small

minFieldSize = 400;
corr_bw_days_mfs = corr_bw_days;
corr_bw_days_mfs(size_place1 < minFieldSize & size_place2 < minFieldSize ) = NaN;

placecell_corr_mean = squeeze(nanmean(corr_bw_days_mfs));
placecell_corr_sem = squeeze(nanstd(corr_bw_days_mfs))./sqrt(numPlaceCells); % in the paper, errorbars were over model iterations

errorbar(smallest_grid,placecell_corr_mean,placecell_corr_sem,placecell_corr_sem,'k')
box off
set(gca,'fontsize',20)
xlabel('grid spacing')
ylabel('correlation')
