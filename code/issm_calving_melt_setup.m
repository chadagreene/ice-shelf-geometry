

load('icemask_composite.mat')
ever_ice = imfill(any(ice,3),8,'holes'); % Take the mean ice mask over the observed period. Taking the sum is faster than mean.  
clear ice year

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
thickness = permute(h5read(fn,'/thickness'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

[X,Y] = meshgrid(x,y); 

ground = ismember(bedmachine_interp('mask',X,Y),[1 2 4]);
bed = bedmachine_interp('bed',X,Y); 
thickness(ground) = bedmachine_interp('thickness',X(ground),Y(ground)); 

thickness(ground & (bed+thickness)<0) = bedmachine_interp('thickness',X(ground & (bed+thickness)<0),Y(ground & (bed+thickness)<0),'method','nearest'); 

thicknessf = filt2(thickness,240,1000,'lp'); 

thickness(ground & (bed+thickness)<0) = thicknessf(ground & (bed+thickness)<0);

ground((thickness+bed)<0) = false; 

sfz = thickness2freeboard(thickness); 
sfz(ground) = bedmachine_interp('surface',X(ground),Y(ground)); 

% vx(~ever_ice) = 32767; 
% vy(~ever_ice) = 32767; 
% thickness(~ever_ice) = 32767; 
% 
vx = int16(vx);
vy = int16(vy); 
thickness = uint16(thickness); 
bed = int16(bed); 
sfz = int16(sfz); 

readme = 'This is the setup for an issm experiment. All velocity and thickness values that are not ever_ice can be set to NaN. Created by issm_calving_melt_setup.mat.';
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','x','y','bed','thickness','sfz','vx','vy','ever_ice','ground','readme','-v7.3')
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup_surface.mat','x','y','sfz','readme','-v7.3') % Only because I neglected to include the surface in the first file, but now issm_calving_melt_setup.mat includes it so no need for a second .mat file. 


