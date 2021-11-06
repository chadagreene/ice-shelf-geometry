
% This script creates the dataset of ice coastline masks, buffered
% inland by varying distances from mean observed coastlines. 
% Run this after icemask_compiler. 
% In 1% increments, this script takes about 2 hours to run. 
% Chad A. Greene, October 2021. 

%% Load data 

load ('icemask_composite.mat')

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 

[X,Y] = meshgrid(x,y); 

ground = ismember(bedmachine_interp('mask',X,Y),[1 2 4]);

%% Calculate max ice mask 

icemax = imfill(any(ice,3),8,'holes'); 

clear ice X Y 
%% Buffer in successive increments of distance from coast

% Range of percent values: 
pct = 0:100; 

% Number of pixels in each ice shelf: 
N = histcounts(mask(icemax & ~ground),.5:181.5);

% Number of pixels to remove in each step: 
Nstep = floor(N/(numel(pct)-1)); 

% Preallocate: 
ice = false(size(icemax,1),size(icemax,2),numel(pct)); 
ice(:,:,1) = icemax; 

for k = 2:length(pct) 
   
   % The most recent ice mask: 
   tmp = ice(:,:,k-1); 
   
   % Distance from coast: 
   dst = bwdist(~tmp); 
   
   for kk = 1:181
      
      % Find which pixels are tidal, correspond to this ice shelf's mask, and were ice in the previous iteration:  
      ind = find(mask==kk & ~ground & tmp); 
      
      % Of the pixels that meet the criteria, sort them by distance to the previous coastline: 
      [~,sort_ind] = sort(dst(ind)); 
   
      % Calve away the N pixels that are closest to the coast:  
      tmp(ind(sort_ind(1:Nstep(kk)))) = false; 
      
   end
   
   % Fix any weird bits: 
   ice(:,:,k) = imfill(tmp,8,'holes'); 
   
   k
end

% Remove any leftovers:  
tmp = ice(:,:,end); 
tmp(~ground) = false; 
ice(:,:,end) = tmp; 

readme = 'maximum ice extents buffered from the coast in 1 percent increments of per-ice-shelf area';  
% save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_buffered.mat','x','y','pct','ice','readme')

%%
% Colors for each year's coastline
col = mat2rgb(year,parula); 

figure
h = imagescn(x,y,ice(:,:,1)); 
bedmachine('gl','color',.5*[1 1 1])
cmocean ice 
hold on
for k = 1:24
   hc(k) = plot(cx{k},cy{k},'color',col(k,:)); 
end
axis off 
set(gcf,'color','k') 
set(gca,'pos',[0 0 1 1])

% Add label on left side 
yr = num2str((1997:2021)');
col = parula(size(yr,1)); 
yr(2:3,:) = []; 
col(2:3,:) = []; 
nc1 = num2str(col(:,1)); 
nc2 = num2str(col(:,2)); 
nc3 = num2str(col(:,3)); 
yr2 = ['\color[rgb]{',nc1(1,:),',',nc2(1,:),',',nc3(1,:),'} ',yr(1,:)];
for k = 2:size(col,1)
   yr2 = [yr2;['\color[rgb]{',nc1(k,:),',',nc2(k,:),',',nc3(k,:),'} ',yr(k,:)]];
end
txt = text(0,0,yr2,'fontsize',7,'units','normalized','vert','bottom','horiz','left'); 

%%

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/coastal_change_maps/icemask_buffer.mp4';

v = VideoWriter(fn,'MPEG-4');
v.FrameRate = 24; 
v.Quality = 100; 
open(v)
   
for k = 1:length(pct)
   h.CData = ice(:,:,k); 
   drawnow
   frame = export_fig('-nocrop','-r600','-painters'); 
   writeVideo(v,frame)
   k
end

close(v)

%% 
function RGB = mat2rgb(val,cmap,limits)
% 
% 
% Chad Greene wrote this, July 2016. 

if nargin==2
   limits = [min(val) max(val)]; 
end

gray = mat2gray(val,limits); 
ind = gray2ind(gray,size(cmap,1));
RGB = cmap(ind+1,:); 

end

