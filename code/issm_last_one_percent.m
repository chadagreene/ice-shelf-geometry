% This script creates a two-panel plot that shows how ISSM responds to the 
% last 1% of ice lost to thinning versus calving. 
% 
% Chad Greene, NASA/JPL January 2022. 


%% Load data

load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_calve_v2.mat') % contains historical then future 
load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_thinningp_v2.mat','thinningp_vel_results') % contains historical then future 
isg = all(ice_mask_results==1,2) & all(ice_vel_results>0,2); 
clear ice_mask_results
%%
%isg = all(ice_vel_results>0,2); 
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

%good = all(ice_vel_results>0,2); 
good = isg; 

%% velocity response 

% From 99% thinned (not 1 m) to 100% calved: 
dvt = ice_vel_results(good,end)-thinningp_vel_results(good,end-1); 
dvc = ice_vel_results(good,end)-ice_vel_results(good,end-1); 

dvpt = 100*dvt./thinningp_vel_results(good,end-1);

dvpc = 100*dvc./ice_vel_results(good,end-1);
% 

calved = ice_vel_results(:,end)==0 & ice_vel_results(:,end-1)>0; 

ax = [-1548391.00    -500346.00     107137.00    1051917.00]; 
ms=0.25; % markersize 

[H,xbm,ybm] = bedmachine_data('thickness',ax(1:2),ax(3:4)); 
mask = bedmachine_data('mask',ax(1:2),ax(3:4)); 
H(mask==0) = nan; 
H = filt2(H,diff(xbm(1:2)),50e3,'lp'); 
H(mask~=3) = nan; 

xl = [-1484087,   -1391952,       -1342045,       -915916,         -789230,      -597280.          ,  -558890,             -562729,        -547373   ,-587999.523965142]; 
yl = [353986,     196587,         131324,          231138,          327113,       296401,           530579,                 757080,        955517.   ,1033734.95933188]; 
tl = {'Evans I.S.';'Carlson Inlet';'Rutford I.S.'; 'Institute I.S.';'Möller I.S.';'Foundation I.S.';{'Support';'Force Gl.'};'Recovery Gl.';'Slessor';'Bailey I.S.'};

xl2 = [-1173129 -685577]; 
yl2 = [737885.43 864572]; 
tl2 = {{'Ronne';'Ice Shelf'};{'Filchner';'Ice Shelf'}}; 

%%

scol = rgb('gold'); 

figure('pos',[ 82.00        451.00        714.00        309.00])

subsubplot(1,2,1,'hpad',0.01) 
hs=fastscatter(x(good),y(good),dvt,'markersize',ms); 
axis tight off
axis(ax)
bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
caxis([-1 1]*500)
cmocean bal
cm = cmocean('bal'); 
ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
uistack(ab,'bottom')
axis(ax)
hold on
[C,hC] = contour(xbm,ybm,H/100,1:1:10,'g');
hC.LineWidth=0.2; 
hC.Color = scol; 
clabel(C,hC,'fontsize',5,'labelspacing',1e5,'color',scol); 
uistack(hC(:),'bottom') 
[C,hC] = contour(xbm,ybm,H/100,11:1:30,'g');
hC.LineWidth=0.2; 
hC.Color = scol; 
% [C,hC] = contour(xbm,ybm,H/100,1:1:30,'g');
% hC.LineWidth=0.2; 
% hC.Color = scol; 
% cbh=clabel(C,hC,'manual'); 
% set(cbh(:),'fontsize',5,'color',scol)
uistack(hC(:),'bottom') 
uistack(hs(:),'top'); 
txt=text(xl,yl,tl,'fontsize',6,'color',0.45*[1 1 1],'fontangle','italic','horiz','center','vert','middle'); 
txt2=text(xl2,yl2,tl2,'fontsize',7,'color',0.45*[1 1 1],'fontangle','italic','horiz','center','vert','middle'); 

cb = colorbar('location','north'); 
set(cb,'xtick',-500:250:500,'fontsize',6)
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
caxis([-1 1]*500)
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
pl.Color = scol; 
uistack(pl,'bottom') 
txt=text(xl,yl,tl,'fontsize',6,'color',0.4*[1 1 1],'fontangle','italic','horiz','center','vert','middle'); 
txt2=text(xl2,yl2,tl2,'fontsize',7,'color',0.4*[1 1 1],'fontangle','italic','horiz','center','vert','middle'); 

ntitle('b','location','nw','fontsize',8,'fontweight','bold')
title('Final 1% calving response','fontsize',8,'fontweight','bold')

%export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_last_one_percent.jpg','-pdf','-r600','-painters','-p0.01')
%%



% 
% 
% % From 99% thinned (not 1 m) to 100% calved: 
% dvt = thinningp_vel_results(good,end-1); 
% 
% dvc = ice_vel_results(good,end-1); 
% 
% calved = ice_vel_results(:,end)==0 & ice_vel_results(:,end-1)>0; 
% 
% ax = [-1548391.00    -500346.00     107137.00    1051917.00]; 
% ms=0.25; % markersize 
% 
% figure('pos',[ 82.00        451.00        714.00        309.00])
% 
% figure
% fastscatter(x(good),y(good),dvt-dvc,'markersize',ms)
% axis tight off
% axis(ax)
% bedmachine('gl','color',rgb('gray'),'linewidth',0.2)
% bedmachine('coast','color',rgb('gray'),'linewidth',0.2)
% caxis([-1 1]*1000)
% cmocean bal
% cm = cmocean('bal'); 
% ab = antbounds('gl','polyshape','edgecolor','none','facecolor',cm(128,:),'facealpha',1);
% uistack(ab,'bottom')
% axis(ax)
% cb = colorbar('location','north'); 
% set(cb,'xtick',-1000:500:1000,'fontsize',6)
% xlabel(cb,'Acceleration (m yr^{-1})','fontsize',6)
% pos = cb.Position; 
% set(cb,'pos',[pos(1)+.02 pos(2)+pos(4)*.8 pos(3)/3 pos(4)/3])
% 
% hold on
% pl = plot(x(calved),y(calved),'.','color',rgb('gray')); 
% pl.MarkerSize=ms;
% pl.Color = 'g'; 
% uistack(pl,'bottom') 

