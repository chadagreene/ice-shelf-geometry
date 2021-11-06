%% Read ITS_LIVE grid

% Read itslive: 
[itsice,x,y] = itslive_data('ice'); 
[X,Y] = meshgrid(x,y); 

yr = [2004 2009 2014]; 
ice = false(size(X,1),size(X,2),3); 

for k = 1:3
   C = load(['moacl',num2str(yr(k)),'.mat']);
   [cx,cy] = ll2ps(C.cllat,C.cllon); 
   ice(:,:,k) = inpolygon_map(X,Y,cx,cy); 
end


%% Visual inspection: 

figure
for k = 1:3 
   subplot(1,3,k) 
   imagescn(x,y,ice(:,:,k)) 
   bedmachine
   ax(k) = gca; 
end
linkaxes(ax,'xy'); 

%% Save it 

readme = 'Ice masks for 2004, 2009, and 2014 from MODIS MOA. Saved by icemask_moa.m'; 
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_moa.mat','x','y','ice','yr','readme','-v7.3')
