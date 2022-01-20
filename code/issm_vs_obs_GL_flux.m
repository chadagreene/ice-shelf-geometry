% This script just creates a quick figure to compare modeled vs
% observed GL flux. 
% Chad Greene, NASA JPL January 2022. 

%%

load('issm_gl_flux_strict.mat')
Dn = load('iceshelves_2008_v2.mat');
C = load('calving_flux_timeseries.mat');

if ~exist('P','var')
   for k = 1:181
      tmpx = Dn.x{k}; 
      tmpy = Dn.y{k}; 
      P(k) = polyshape(tmpx,tmpy); 
   end
end
Asc = sqrt(Dn.area_km2); 
Asc = 90*Asc/max(Asc); 

[~,lon] = ps2ll(Dn.x_center,Dn.y_center);
col = mat2rgb(lon,cmocean('phase'),[-1 1]*180); 
col(183,:) = [0 0 0];

%%

Dn.name(182)={'Other'}; 
Dn.name(183)={'Antarctica'};  
figure('pos',[40.00        534.00        316.00        341.00])

scatter(C.GL_ss(1:181),glf_0(1:181),Asc,lon,'filled')
cmocean phase
caxis([-1 1]*180)
set(gca,'xscale','log','yscale','log')
axis tight
daspect([1 1 1]) 
axis([1e-1 250 1e-1 250]) 

hold on
pl=plot(xlim,ylim,'color',rgb('gray'),'linewidth',0.25); 
uistack(pl,'bottom') 

xlabel('Measured GL flux (Gt yr^{-1})','fontsize',7)
ylabel('Modeled control GL flux (Gt yr^{-1})','fontsize',7)

txt = text(C.GL_ss(1:181),glf_0(1:181),Dn.name(1:181),'color',.4*[1 1 1],'fontangle','italic','horiz','center','vert','bot','fontsize',5,'clipping','on'); 

set(gca,'fontsize',7)

gp = plotboxpos(gca);

axes('position',[gp(1)+(.7)*gp(3) gp(2) .3*gp(3) .3*gp(4)])
 
 
hold on
for k = 1:181
   plot(P(k),'facecolor',col(k,:),'facealpha',1,'edgecolor','none')
end

hant = antbounds('coast','polyshape','facecolor','w','facealpha',1,'edgecolor','none');
uistack(hant,'bottom'); 
axis tight off
bedmachine('coast','color',0.3*[1 1 1],'linewidth',0.1)
bedmachine('gl','color',0.3*[1 1 1],'linewidth',0.1)
gp2 = plotboxpos(gca); 
set(gca,'xcolor','none','ycolor','none','pos',[gp(1)+(.7)*gp(3) gp(2) gp2(3) gp2(4)])
set(gca,'fontsize',7)

% export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_vs_obs_GL_flux.jpg','-pdf','-r600','-painters')

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
