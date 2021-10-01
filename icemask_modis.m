% This script reads and interpolates Alex Fraser's annual modis coastline
% picks to the ITS_LIVE 240 m grid. 
% 
% Chad Greene, September 2021. 

%% Load data

% Read itslive grid: 
[itsice,x,y] = itslive_data('ice'); 
[X,Y] = meshgrid(x,y); 

% The public data: 
[fi,xfi,yfi,years] = fastice_data_autumn;

tmp = flipud(permute(ncread('coast_2021_old_resample_qc.nc','coastline'),[2 1])); 

fi = cat(3,fi,tmp); 
years = [years;2021]; 

%% Interpolate to the ITS_LIVE grid: 

ice = false(length(y),length(x),length(years)); 

for k = 1:length(years)
   % fi values of 1 2 3 correspond to continent, islands, ice shelf. 
   ice(:,:,k) = interp2(xfi,yfi,ismember(fi(:,:,k),1:3),X,Y,'nearest'); 
   
   % fill any holes: 
   ice(:,:,k) = imfill(ice(:,:,k),'holes'); 
end

readme = 'ice mask from Alex Frasers digitization of modis data, interpolated to the its_live grid by icemask_modis.m'; 

% save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_modis.mat','x','y','years','ice','readme','-v7.3')
