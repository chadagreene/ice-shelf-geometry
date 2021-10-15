

load('icemask_composite.mat')
ever_ice = imfill(any(ice,3),8,'holes'); % Take the mean ice mask over the observed period. Taking the sum is faster than mean.  
clear ice year

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-15.h5'; 
thickness = permute(h5read(fn,'/thickness'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

[X,Y] = meshgrid(x,y); 

ground = ismember(bedmachine_interp('mask',X,Y),[1 2 4]);
bed = bedmachine_interp('bed',X,Y); 

vx(~ever_ice) = 32767; 
vy(~ever_ice) = 32767; 
thickness(~ever_ice) = 32767; 

vx = int16(vx);
vy = int16(vy); 
thickness = uint16(thickness); 
bed = int16(bed); 

readme = 'This is the setup for an issm experiment. Created by issm_calving_melt_setup.mat.';
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','x','y','bed','thickness','vx','vy','ever_ice','ground','readme','-v7.3')
