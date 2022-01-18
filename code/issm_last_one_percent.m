


load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_calve_v2.mat','x','y','ice_vel_results') % contains historical then future 
load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_thinningp_v2.mat','thinningp_vel_results') % contains historical then future 

%%
isg = all(ice_vel_results>0,2); 
% For calving experiment (past) assume no acceleration occurred on grounded ice since last timestep: 
% And start over for hypothetical future experiment, but don't let the past affect the future: 
for k = [3:25 27:126]
   
   ind = (isg & ice_vel_results(:,k)==0) | ice_vel_results(:,k)>10e3; 
   ice_vel_results(ind,k) = ice_vel_results(ind,k-1); 
   
end

% Do the same for future thinning
for k = 2:100
   
   ind = (isg & thinningp_vel_results(:,k)==0) | thinningp_vel_results(:,k)>10e3; 
   thinningp_vel_results(ind,k) = thinningp_vel_results(ind,k-1); 
   
end

good = all(ice_vel_results>0,2); 

%%

% From 99% thinned (not 1 m) to 100% calved: 
dv = ice_vel_results(good,end)-thinningp_vel_results(good,end-1); 
dvpt = 100*dv./thinningp_vel_results(good,end-1);

dv = ice_vel_results(good,end)-ice_vel_results(good,end-1); 
dvpc = 100*dv./ice_vel_results(good,end-1);

calved = ice_vel_results(:,end)==0 & ice_vel_results(:,end-1)>0; 

ax = [-1548391.00    -500346.00     107137.00    1051917.00]; 
ms=0.25; % markersize 

figure('pos',[ 82.00        451.00        714.00        309.00])

subsubplot(1,2,1,'hpad',0.01) 
fastscatter(x(good),y(good),dvpt,'markersize',ms)
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*100)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
cb = colorbar('location','north'); 
set(cb,'xtick',-100:50:100,'fontsize',6)
xlabel(cb,'Acceleration (%)','fontsize',6)
pos = cb.Position; 
set(cb,'pos',[pos(1)+.02 pos(2)+pos(4)*.8 pos(3)/3 pos(4)/3])

ntitle('a','location','nw','fontsize',8,'fontweight','bold')
title('Final 1% thinning response','fontsize',8,'fontweight','bold')

subsubplot(1,2,2,'hpad',0.01) 
fastscatter(x(good),y(good),dvpc,'markersize',ms)
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*100)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
[hs1,hs2]=scalebarps('location','se','fontsize',6,'length',100); 
hs2.LineWidth = 1; 
hold on
pl = plot(x(calved),y(calved),'.','color',rgb('gray')); 
pl.MarkerSize=ms;
pl.Color = 'g'; 
uistack(pl,'bottom') 

ntitle('b','location','nw','fontsize',8,'fontweight','bold')
title('Final 1% calving response','fontsize',8,'fontweight','bold')

%% velocity response 


% From 99% thinned (not 1 m) to 100% calved: 
dvt = ice_vel_results(good,end)-thinningp_vel_results(good,end-1); 

dvc = ice_vel_results(good,end)-ice_vel_results(good,end-1); 

calved = ice_vel_results(:,end)==0 & ice_vel_results(:,end-1)>0; 

ax = [-1548391.00    -500346.00     107137.00    1051917.00]; 
ms=0.25; % markersize 

figure('pos',[ 82.00        451.00        714.00        309.00])

subsubplot(1,2,1,'hpad',0.01) 
fastscatter(x(good),y(good),dvt,'markersize',ms)
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*1000)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
cb = colorbar('location','north'); 
set(cb,'xtick',-1000:500:1000,'fontsize',6)
xlabel(cb,'Acceleration (m yr^{-1})','fontsize',6)
pos = cb.Position; 
set(cb,'pos',[pos(1)+.02 pos(2)+pos(4)*.8 pos(3)/3 pos(4)/3])

ntitle('a','location','nw','fontsize',8,'fontweight','bold')
title('Final 1% thinning response','fontsize',8,'fontweight','bold')

subsubplot(1,2,2,'hpad',0.01) 
fastscatter(x(good),y(good),dvc,'markersize',ms)
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*1000)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
[hs1,hs2]=scalebarps('location','se','fontsize',6,'length',100); 
hs2.LineWidth = 1; 
hold on
pl = plot(x(calved),y(calved),'.','color',rgb('gray')); 
pl.MarkerSize=ms;
pl.Color = 'g'; 
uistack(pl,'bottom') 

ntitle('b','location','nw','fontsize',8,'fontweight','bold')
title('Final 1% calving response','fontsize',8,'fontweight','bold')

%%





% From 99% thinned (not 1 m) to 100% calved: 
dvt = thinningp_vel_results(good,end-1); 

dvc = ice_vel_results(good,end-1); 

calved = ice_vel_results(:,end)==0 & ice_vel_results(:,end-1)>0; 

ax = [-1548391.00    -500346.00     107137.00    1051917.00]; 
ms=0.25; % markersize 

figure('pos',[ 82.00        451.00        714.00        309.00])

figure
fastscatter(x(good),y(good),dvt-dvc,'markersize',ms)
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*1000)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
cb = colorbar('location','north'); 
set(cb,'xtick',-1000:500:1000,'fontsize',6)
xlabel(cb,'Acceleration (m yr^{-1})','fontsize',6)
pos = cb.Position; 
set(cb,'pos',[pos(1)+.02 pos(2)+pos(4)*.8 pos(3)/3 pos(4)/3])

hold on
pl = plot(x(calved),y(calved),'.','color',rgb('gray')); 
pl.MarkerSize=ms;
pl.Color = 'g'; 
uistack(pl,'bottom') 

