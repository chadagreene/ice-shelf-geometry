% This script reads ice masks from Sentinel 1 a. 
% The masks are on a 1 km grid and correspond to days 56-80 of each year. 
% The masks were generated by Alex Fraser from the data at seaice.dk. 
% Alex Fraser saved them as .xcf (GIMP) files, then I converted them to 
% tif files by the same name. 
%
% This script simply reads the masks and interpolates to the ITS_LIVE 
% velocity grid. 
% 
% Chad Greene, Sept 2021. 

%% Read ITS_LIVE grid
% (And convert to spherical ps70)

% Read itslive: 
[itsice,x,y] = itslive_data('ice'); 

[X,Y] = meshgrid(x,y); 
[Lat,Lon] = ps2ll(X,Y); 

% ITS_LIVE grid locations on the spherical ps70 grid: 
[X70,Y70] = ll2ps(Lat,Lon,...
   'EarthRadius',6370997.0,...
   'Eccentricity',0,...
   'TrueLat',-70); 

%% Read annual ice masks: 

yr = 2015:2021; % years

% Read first image: 
tmp = imread('2015_056-080.tif'); 

% Preallocate: 
ices = zeros(size(tmp,1),size(tmp,2),length(yr)); 

% Populate the first year with data from the first image: 
ices(:,:,1) = tmp(:,:,1)>0; 

% Loop through the remaining images: 
for k = 2:length(yr) 
   tmp = imread([num2str(yr(k)),'_056-080.tif']); 
   ices(:,:,k) = tmp(:,:,1)>0; 
end

%%

% Grid cell center of upper left corner: (according to email with Roberto Saldo) 
[xs0,ys0] = ll2ps(-39.219500, -42.240700,...
   'EarthRadius',6370997.0,...
   'Eccentricity',0,...
   'TrueLat',-70); 

% The grid is equally spaced 1 km in spherical ps70 coordinates: 
xs70 = xs0 + 1000*(1:size(ices,2)); 
ys70 = ys0 - 1000*(1:size(ices,1));

% Preallocate ice mask on itslive grid: 
ice = false(size(Lat,1),size(Lat,2),length(yr)); 

% Interpolate spherical ps70 to itslive grid: 
for k = 1:length(yr) 
   ice(:,:,k) = interp2(xs70,ys70,ices(:,:,k),X70,Y70,'nearest'); 
   
   % fill any holes: 
   ice(:,:,k) = imfill(ice(:,:,k),'holes'); 
end

readme = 'written by icemask_s1a.m. Ice mask from Alex Frasers digitization of the Sentinel 1a data on a 1 km grid.';
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_s1a.mat','x','y','ice','readme','yr','-v7.3')

%% Verify alignment 
% The spherical ps70 is strange, but I think it lines up. 

[clat,clon] = antbounds_data('coast'); 
[cx70,cy70] = ll2ps(clat,clon,...
   'EarthRadius',6370997.0,...
   'Eccentricity',0,...
   'TrueLat',-70);

figure
imagescn(xs70,ys70,tmp(:,:,1))
hold on
plot(cx70,cy70,'k')
daspect([1 1 1])
