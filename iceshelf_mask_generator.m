
% This script creates a mask of Mouginot's iceshelves_2008_v2 on 
% the 240 m ITS_LIVE grid. It dilates the extents of all ice shelves 
% to allow for changing ice shelf geometry. 
% 
% This script saves iceshelf_mask.mat, which contains 

% Chad A. Greene NASA/JPL August 2021. 

%% Requirements: 

% * inpolygon_map.m function by Chad Greene
% * iceshelves_2008_v2.mat from Antarctic Boundaries (on file exchange) 
% * itslive_data function (on github) 
% * itslive mosaic data 
% * image processing toolbox 
% * mapping toolbox (for inpolygon_map function) 

%% Load data 

% Just the itslive grid: 
[~,x,y] = itslive_data('ocean'); 
[X,Y] = meshgrid(x,y); 

% Mouginot's ice shelf mask: 
D = load('iceshelves_2008_v2.mat'); % requires antbounds for AMT. 
names = D.name; 

%% Mask each ice shelf 
% This will take a few minutes. 

iceshelves = zeros(size(X),'uint8'); 

f = waitbar(0);
for k = 1:length(names) 
   iceshelves(inpolygon_map(X,Y,D.x{k},D.y{k})) = k; 
   waitbar(k/181,f)
end
close(f)

figure
imagescn(x,y,iceshelves)
bedmachine

%% Dilate the ice shelves

% Start with the initial data: 
iceshelves_dilated = iceshelves; 

% Dilation element: 
st = strel('disk',4); % Each dilation is by 4 px ~ 1 km 

% Do the dilation in multiple steps so ice shelves with a large numbered
% index (near the end of the alphabet) don't overtake smaller index numbers
% (because the imdilate function prioritizes the largest number).
for k = 1:104 % Go to 104*4*0.240 ~ 100 km
   
   % Make a temporary record of the current mask: 
   iceshelves_dilated_old = iceshelves_dilated;
   
   % Dilate all the ice shelves: 
   iceshelves_dilated = imdilate(iceshelves_dilated,st); 

   % If any ice shelf was present in the old mask, make sure the new mask
   % didn't overwrite it: 
   iceshelves_dilated(iceshelves_dilated_old>0) = iceshelves_dilated_old(iceshelves_dilated_old>0);

   k % just a counter 
end

figure
imagescn(x,y,iceshelves_dilated)
bedmachine

readme = 'Mouginot et als 2008 v2 ice shelf outlines applied to the ITS_LIVE grid, written by iceshelf_mask_generator.m. Ice shelves are dilated by 100 km.'; 

% save('/Users/cgreene/Documents/data/my_data/iceshelf_mask.mat','x','y','iceshelves','iceshelves_dilated','names','readme'); 

