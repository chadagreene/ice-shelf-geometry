% This script provides a quick comparison to figure S12 Gudmudsson's 2019
% paper. Our results are very similar, although it may not look like it at
% first glance becuase Gudmundsson's colormap is a bit more sensitive to
% tiny variations near the zero value. 
% 
% Chad Greene, NASA JPL January 2022. 

%% load data

load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_thickness_v2.mat','x','y','thickness_vel_results') % contains historical then future 

good = all(thickness_vel_results>0,2); 
good(isopenocean(x,y)) = false; 

%% Plot up the velocity response
% Compare 2017 to 2010 

dv = thickness_vel_results(good,26)-thickness_vel_results(good,19); 
dvp = 100*dv./thickness_vel_results(good,19);

figure('pos',[82.00        457.00        357.00        286.00])
fastscatter(x(good),y(good),dvp,'markersize',1)
axis tight off
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*10)
cmocean bal
cb = colorbarps('fontsize',6); 
set(cb,'xtick',-10:5:10)
xlabel(cb,'Change in speed (%)')

cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')

% export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_gudmundsson_compare.jpg','-r600','-painters')
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/fig_ED7.eps','contenttype','vector')


%%

good = good & isiceshelf(x,y); 

dv = thickness_vel_results(good,26)-thickness_vel_results(good,19); 
dvp = 100*dv./thickness_vel_results(good,19);

figure('pos',[82.00        457.00        357.00        286.00])
fastscatter(x(good),y(good),dvp,'markersize',1)
axis tight off
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*10)
cmocean bal
cb = colorbarps('fontsize',6); 
set(cb,'xtick',-10:5:10)
xlabel(cb,'Change in speed (%)')

cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')

cm = [flipud(brewermap(128,'blues'));
brewermap(128,'YlOrRd')]; 

colormap(cm)

% export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_gudmundsson_compare_brewermap.jpg','-r600','-painters')

