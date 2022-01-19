% This script creates an animation of Antarctica's gl flux response 
% to ice shelf collapse. 
% Run this after issm_thickness_response_analysis.m
% 
% Chad Greene, NASA/JPL January 2022. 

%% Load data 

%D = load('/Users/cgreene/Documents/MATLAB/DEM_generation/iceshelf_thickness_cube_1992-2018.mat','x','y','year','ground'); 
load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_calve_v2.mat','ice_mask_results') % contains control, then historical, then future

G = load('issm_gl_flux_strict.mat'); % from issm_thickness_response_analysis.m
%%

isg = all(ice_mask_results==1,2); 

clear ice_mask_results

% For calving experiment (past) assume no acceleration occurred on grounded ice since last timestep: 
% And start over for hypothetical future experiment, but don't let the past affect the future: 
for k = [3:25 27:126]
   
   ind = (isg & ice_vel_results(:,k)==0) | ice_vel_results(:,k)>10e3; 
   ice_vel_results(ind,k) = ice_vel_results(ind,k-1); 

end

%% Initialize the figure:

cm = cmocean('bal'); 

k=1;

figure('color',rgb('gray'))
tmp = ice_vel_results(:,25+k)-ice_vel_results(:,26); 
good = ice_vel_results(:,25+k)>0; 
h = fastscatter(x,y,tmp); 
hold on
box off 
axis tight off
bedmachine('gl','color',rgb('gray'),'linewidth',0.3)
set(gca,'pos',[0 0 1 1])
cmocean bal
caxis([-1 1]*1000)
cb = colorbarps; 
xlabel(cb,'Change in Ice Speed (m/yr)','fontsize',7)
set(cb,'fontsize',7)
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
drawnow
axmap=gca; 

ax = axes; 
pl1 = plot(0:100,G.glf_cf(end,:),'color',rgb('dark red'),'linewidth',1); 
hold on
pl2 = plot(100,G.glf_cf(end,end),'o','markersize',4,'markerfacecolor',rgb('dark red'),'color',rgb('dark red'),'linewidth',1); 
set(ax,'pos',[.45 .45 .32 .38],'color','none','fontsize',7,'xtick',0:10:100)
xlabel('percent ice shelf collapse','fontsize',7)
box off
axis([0 100 2100 5200])
title 'Grounding Line Flux (Gt/yr)'

axes(axmap) % makes map axes current 
uistack(axmap,'bottom') 

%% Write the video

vd = VideoWriter(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/animations/issm_hypothetical_calving_response.mp4'],'MPEG-4'); 
vd.FrameRate=16; %
vd.Quality=100; 
open(vd)

for k=1:101
   delete(h); 
   tmp = ice_vel_results(:,25+k)-ice_vel_results(:,26); 
   good = ice_vel_results(:,25+k)>0; 
   h = fastscatter(x(good),y(good),tmp(good)); 
   uistack(h,'bottom')
   uistack(ab,'bottom'); 

   set(pl1,'XData',0:(k-1),'YData',G.glf_cf(end,1:k))
   set(pl2,'XData',k-1,'YData',G.glf_cf(end,k))
   drawnow
   fr = export_fig('-nocrop','-r600','-painters'); % can use getframe(gcf) instead, but export_fig offers more control over resolution, etc.  
   writeVideo(vd,fr)
   k
end
close(vd)

clear good cb fr h k tmp vd
