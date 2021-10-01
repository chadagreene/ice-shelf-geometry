%% Read ITS_LIVE grid

% Read itslive: 
[itsice,x,y] = itslive_data('ice'); 
[X,Y] = meshgrid(x,y); 

% The old way: ( Slower and less accurate): 
% [Iramp,xramp,yramp] = geoimread('ramp2_dem_osu91a200m.tif'); % The osu91a version of the file is referenced to the geoid whereas the RAM2_DEM.tif is referenced to wgs84. Use OSU.  
% ice = interp2(xramp,yramp,Iramp>0,X,Y,'nearest'); 

% Read shapefiles: 
S1997 = shaperead('after_coast_continuous.shp'); 
S2000 = shaperead('cst2000line.shp');

% Remove islands from the 2000 sata: 
S2000 = S2000([S2000.RPOLY_]==1 & [S2000.LPOLY_]==2);

% Create a continuous polygon: 
[tx,ty] = polymerge([S2000.X],[S2000.Y]); 

% Fix one out-of-order segment: 
tx = [tx(32725:-1:1),tx(32726:end)]; 
ty = [ty(32725:-1:1),ty(32726:end)]; 

isf = isfinite(tx); 

% Create masks: 
ice1997 = inpolygon_map(X,Y,S1997.X(1:end-1),S1997.Y(1:end-1)); 
ice2000 = inpolygon_map(X,Y,tx(isf),ty(isf)); 

% Fill Holes: 
ice1997 = imfill(ice1997,'holes'); 
ice2000 = imfill(ice2000,'holes'); 

%% Visual inspection: 

figure
subplot(1,2,1) 
imagescn(x,y,ice1997) 
bedmachine
hold on
plot(S1997.X,S1997.Y,'r')
ax(1) = gca; 

subplot(1,2,2) 
imagescn(x,y,ice2000) 
bedmachine
hold on
plot(tx(isf),ty(isf),'r')
ax(2) = gca; 

linkaxes(ax,'xy')

%% Save it 

ice = cat(3,ice1997,ice2000); 
yr = [1997 2000];
readme = 'Ice mask for 1997 and 2000 from after_coast_continuous.shp and cst2000line.shp, respectively. Saved by icemask_ramp2.m'; 

% save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_ramp2b.mat','x','y','ice','yr','readme','-v7.3')
