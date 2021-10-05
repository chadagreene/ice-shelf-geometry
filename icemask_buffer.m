
% This script creates the dataset of ice coastline masks, buffered
% inland by varying distances from mean observed coastlines. 
% Run this after icemask_compiler. 

%% Load data 

load ('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_composite.mat')


[X,Y] = meshgrid(x,y); 

ground = ismember(bedmachine_interp('mask',X,Y),[1 2 4]);

%% Calculate mean ice mask 
% The mean function struggles on this big of a dataset if it's double, so 
% sum it up as single and then divide by 24. 

icem = logical(sum(ice,3)/24);

icem = imfill(icem,8,'holes'); 

figure
imagescn(x,y,icem) 
bedmachine('gl','color',[.5 .5 .5])

%% Buffer 

d = [0:.25:10 11:20];

