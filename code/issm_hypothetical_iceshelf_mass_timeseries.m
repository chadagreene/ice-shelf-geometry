% This script creates ice shelf mass (time) series for hypothetical
% situations where 
% 1. ice shelves continue their recent trends for the next 50 years, and 
% 2. ice shelves instantaneously lose a given percentage of their current area. 
% 
% This script takes about 5 or 10 minutes to run on my laptop. 
% Chad A. Greene, NASA Jet Propulsion Laboratory, January 2022. 

%% Load data 

load('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_buffered.mat')
load('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','thickness','ground','ever_ice')
mask = permute(ncread('/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5','iceshelf_mask'),[2 1]); 

%%

[X,Y] = meshgrid(x,y); 

% Determine the mass of ice in each grid cell: 
[Lat,~] = ps2ll(X,Y); 
A = (diff(x(1:2)) ./ psdistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  
F = 917.*A.*1e-12; % a multiplication factor
IceMass = 917.*double(thickness).*A.*1e-12; % ice mass in each grid cell (Gt)
clear Lat A thickness X Y

%%

ice2 = cube2rect(ice,ever_ice & ~ground); 
clear ice 

IceMass2 = cube2rect(IceMass,ever_ice & ~ground); 
mask2 = cube2rect(mask,ever_ice & ~ground); 
M2 = double(ice2).*IceMass2;

%%

M_calving_future = nan(101,181); 
for k = 1:181
   M_calving_future(:,k) = sum(M2(:,mask2==k),2); 
end
M_calving_future(:,182) = sum(M2(:,mask2==0),2); 
M_calving_future(:,183) = sum(M2,2); 

%%

clear ice2 M2 IceMass* mask2
load('/Users/cgreene/Documents/MATLAB/DEM_generation/iceshelf_thickness_cube_future.mat')

H2 = cube2rect(H,always_ice & ~ground); 
F2 = cube2rect(F,always_ice & ~ground); 
mask2 = cube2rect(mask,always_ice & ~ground); 
clear F H mask

M2 = double(H2).*F2; 

M_thinning_future = nan(26,181); 
for k = 1:181
   M_thinning_future(:,k) = sum(M2(:,mask2==k),2); 
end
M_thinning_future(:,182) = sum(M2(:,mask2==0),2); 
M_thinning_future(:,183) = sum(M2,2); 

%%

readme = 'mass time series in Gt for ice shelves. Written by issm_hypothetical_iceshelf_mass_timeseries.m'; 
save('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/hypothetical_iceshelf_mass.mat','pct','M_calving_future','future_years','M_thinning_future','readme')

