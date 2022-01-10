% This script converts the output of ISSM runs into grounding line flux
% for each ice shelf. 
% 
% 
% Chad Greene, NASA/JPL January 2022. 

%% Load data 

load('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','ground','thickness','vx','vy')
D = load('/Users/cgreene/Documents/MATLAB/DEM_generation/iceshelf_thickness_cube_1992-2018.mat','x','y','year'); 
I = load('icemask_composite.mat','year'); 
load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_calve_v2.mat') % contains control, then historical, then future
load('/Users/cgreene/Documents/MATLAB/DEM_generation/ice_vel_results_thickness_v2.mat') % contains historical then future 

isg = all(ice_mask_results==1,2); 

H0 = thickness_H_results(:,1); % only take the first thickness bc that's the control. 
clear thickness_H_results ice_mask_results

%% Constrain Island velocities after calving

% For calving experiment (past) assume no acceleration occurred on grounded ice since last timestep: 
% And start over for hypothetical future experiment, but don't let the past affect the future: 
for k = [3:25 27:126]
   
   ind = (isg & ice_vel_results(:,k)==0) | ice_vel_results(:,k)>10e3; 
   ice_vel_results(ind,k) = ice_vel_results(ind,k-1); 
   
end

%% Trim datasets 
% To ensure we're fitting a surface to the same nodes in every experiment, 
% we now remove any nodes that aren't always there. 

keep = any(ice_vel_results>0,2) & all(thickness_vel_results>0,2); 
ice_vel_results = ice_vel_results(keep,:); 
thickness_vel_results = thickness_vel_results(keep,:); 
x = x(keep); 
y = y(keep); 
H0 = H0(keep); 

clear keep 

%% Separate results into a matrix for each experiment

v0 = ice_vel_results(:,1); % control run
vcp = ice_vel_results(:,2:25); % calving past
vcf = ice_vel_results(:,26:end); % calving future

vtp = thickness_vel_results(:,1:26); % thickness past
vtf = thickness_vel_results(:,27:end); % thickness future
clear thickness_* ice_vel* ind isg

%% Define the grounding line: 

% Our original input velocity and thickness at each node: 
v = interp2(D.x,D.y,hypot(double(vx),double(vy)),x,y); 
H = interp2(D.x,D.y,double(thickness),x,y); 

% Convert input ground mask to an outline: 
[xgl,ygl] = mask2outline(D.x,D.y,ground); 

% Trim the grounding line to only include islands that have ISSM nodes inside them:  
[xglc,yglc] = polysplit(xgl,ygl); 
good = false(size(xglc)); 
for k=1:length(xglc)
   good(k) = any(fast_inpolygon(x,y,xglc{k},yglc{k}));
end
[xgl,ygl] = polyjoin(xglc(good),yglc(good));

% Calculate thickness and velocity along the GL from gridded data: 
Hgl = interp2(D.x,D.y,double(thickness),xgl,ygl); 
vxgl = interp2(D.x,D.y,double(vx),xgl,ygl); 
vygl = interp2(D.x,D.y,double(vy),xgl,ygl); 

% Get unit vectors of velocity: 
vgl = hypot(vxgl,vygl);
vxgl_unit = vxgl./vgl;  
vygl_unit = vygl./vgl; 

clear vx vy thickness good

% For reference calculate GL flux using the original input data:  
flux = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*vgl.*Hgl*917*1e-12;
sum(flux,'omitnan')

%% Get GL thickness and velocity from ISSM: 

% Set up an interpolant that is the difference between ISSM thickness and input
% gridded thickness. The rationale here is that fitting a surface to the difference 
% means smaller range and hopefully smaller interpolation errors: 
keep = v0>0; 
F = scatteredInterpolant(x(keep),y(keep),H0(keep)-H(keep),'natural'); 

% Get the reference thickness along the GL from ISSM: 
H0gl = F(xgl,ygl)+Hgl; 

% Update the interpolant for reference GL velocity (faster than creating a new interpolant from scratch)  
F.Values = v0(keep)-v(keep); 
v0gl = F(xgl,ygl)+vgl; 

% Calculate GL flux for the reference run: 
flux0 = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*v0gl.*H0gl*917*1e-12; 
sum(flux0,'omitnan')

%% Calving Past

keep = all(vcp>0,2); 

fluxcp = nan(numel(flux0),size(vcp,2)); 

for k=1:size(fluxcp,2)
   
   if k==1
      F = scatteredInterpolant(x(keep),y(keep),vcp(keep,k)-v0(keep),'natural'); 
   else
      % Update the interpolant: 
      F.Values = vcp(keep,k)-v0(keep); 
   end
   
   % Calculate velocity along the gl: 
   vtmp = F(xgl,ygl) + v0gl; 

   % Get flux at every grid point: 
   fluxcp(:,k) = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*(vtmp).*H0gl*917*1e-12; 

   k
end

figure
plot(I.year,-sum(fluxcp,'omitnan'))
box off 
axis tight
%hline(-sum(flux0,'omitnan'),'color',rgb('gray'))
title 'calving response' 
ylabel 'gl flux (Gt/yr)'

%% Calving future: 

keep = all(vcf>0,2); 

fluxcf = nan(numel(flux0),size(vcf,2)); 

for k=1:101
   
   if k==1
      F = scatteredInterpolant(x(keep),y(keep),vcf(keep,k)-v0(keep),'natural'); 
   else
      % Update the interpolant: 
      F.Values = vcf(keep,k)-v0(keep); 
   end
   
   % Calculate velocity along the gl: 
   vtmp = F(xgl,ygl) + v0gl; 

   % Get flux at every grid point: 
   fluxcf(:,k) = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*vtmp.*H0gl*917*1e-12; 

   k
end

figure
plot(0:100,-sum(fluxcf,'omitnan'))
hold on
box off 
axis tight
ylabel 'gl flux (Gt/yr)'
xlabel 'percent ice shelf loss'
title 'pan-Antarctic response to ice shelf collapse'
% export_fig issm_response_calving_hypothetical.png -r500 -p0.01

%% Thinning past

keep = all(vtp>0,2); 
fluxtp = nan(numel(flux0),size(vtp,2)); 

for k=1:size(fluxtp,2)
   
   if k==1
      F = scatteredInterpolant(x(keep),y(keep),vtp(keep,k)-v0(keep),'natural'); 
   else
      % Update the interpolant: 
      F.Values = vtp(keep,k)-v0(keep); 
   end
   
   % Calculate velocity along the gl: 
   vtmp = F(xgl,ygl)+v0gl; 

   % Get flux at every grid point: 
   fluxtp(:,k) = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*vtmp.*H0gl*917*1e-12; 

   k
end

figure
plot(D.year,-sum(fluxtp,'omitnan')-mean(-sum(fluxtp,'omitnan') ))
box off 
axis tight
hold on
pl = plot(I.year,-sum(fluxcp,'omitnan') - mean(-sum(fluxcp,'omitnan') ));
legend('thickness','calving','location','northwest')
ylabel 'GL flux anomaly (Gt/yr)'

% export_fig issm_response_thinning_calving_past.png -r500 -p0.01
%% Thinning future

load('iceshelf_thickness_cube_future.mat','future_years') 

keep = all(vtf>0,2); 
fluxtf = nan(numel(flux0),size(vtf,2)); 

for k=1:size(fluxtf,2)
   
   if k==1
      F = scatteredInterpolant(x(keep),y(keep),vtf(keep,k)-v0(keep),'natural'); 
   else
      % Update the interpolant: 
      F.Values = vtf(keep,k)-v0(keep); 
   end
   
   % Calculate velocity along the gl: 
   vtmp = F(xgl,ygl) + v0gl; 

   % Get flux at every grid point: 
   fluxtf(:,k) = (gradient(xgl).*vygl_unit - gradient(ygl).*vxgl_unit).*vtmp.*H0gl*917*1e-12; 

   k
end

figure
plot(future_years,-sum(fluxtf,'omitnan'))
box off 
axis tight

%%
Dn = load('iceshelves_2008_v2.mat');

fn = 'extruded_antarctica_2021-10-18.h5';
shelf = interp2(ncread(fn,'x'),ncread(fn,'y'),permute(ncread(fn,'iceshelf_mask'),[2 1]),xgl,ygl,'nearest'); 

glf_cp = nan(183,size(fluxcp,2)); 
for k = 1:181
   glf_cp(k,:) = -sum(fluxcp(shelf==k,:),1,'omitnan'); 
end
glf_cp(182,:) = -sum(fluxcp(~ismember(shelf,1:181),:),1,'omitnan'); 
glf_cp(183,:) = -sum(fluxcp,1,'omitnan'); 


glf_tp = nan(183,size(fluxtp,2)); 
for k = 1:181
   glf_tp(k,:) = -sum(fluxtp(shelf==k,:),1,'omitnan'); 
end
glf_tp(182,:) = -sum(fluxtp(~ismember(shelf,1:181),:),1,'omitnan'); 
glf_tp(183,:) = -sum(fluxtp,1,'omitnan'); 


glf_cf = nan(183,size(fluxcf,2)); 
for k = 1:181
   glf_cf(k,:) = -sum(fluxcf(shelf==k,:),1,'omitnan'); 
end
glf_cf(182,:) = -sum(fluxcf(~ismember(shelf,1:181),:),1,'omitnan'); 
glf_cf(183,:) = -sum(fluxcf,1,'omitnan'); 


glf_tf = nan(183,size(fluxtf,2)); 
for k = 1:181
   glf_tf(k,:) = -sum(fluxtf(shelf==k,:),1,'omitnan'); 
end
glf_tf(182,:) = -sum(fluxtf(~ismember(shelf,1:181),:),1,'omitnan'); 
glf_tf(183,:) = -sum(fluxtf,1,'omitnan'); 

%% Individual response to past calving

% Line colors: 
col = double(lab2rgb([80*rand(183,1)+5 128*2*(rand(183,1)-0.5) 128*2*(rand(183,1)-0.5)],'OutputType','uint8'))/255;
col(end,:) = [0 0 0];


figure('pos',[39.00        562.00        603.00        536.00]) 
% Individual response to past thinning 
subsubplot(1,2,1,'hpad',0.05) 
hold on
for k=1:181
   pl(k) = plot(D.year,glf_tp(k,:)-glf_tp(k,1),'color',col(k,:)); 
   if k<=181
      text(D.year(end),glf_tp(k,end)-glf_tp(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:))
   end
end
pl(183) = plot(D.year,glf_tp(183,:)-glf_tp(183,1),'color',col(183,:)); 
pl(183).LineWidth = 1; 
text(D.year(end),glf_tp(183,end)-glf_tp(183,1),'Antarctica','fontsize',6,'horiz','left','vert','middle','color',col(183,:),'fontweight','bold')
box off
axis tight
ylabel 'change in GL flux (Gt/yr)'
ntitle('Response to Thickness Change','fontsize',8)
ylim([-2.5 19])
set(gca,'fontsize',7)

subsubplot(1,2,2,'hpad',0.05) 
hold on
for k=1:181
   pl(k) = plot(I.year,glf_cp(k,:)-glf_cp(k,1),'color',col(k,:)); 
   if k<=181
      text(I.year(end),glf_cp(k,end)-glf_cp(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:))
   end
end
pl(183) = plot(I.year,glf_cp(183,:)-glf_cp(183,1),'color',col(183,:)); 
pl(183).LineWidth = 1; 
text(I.year(end),glf_cp(183,end)-glf_cp(183,1),'Antarctica','fontsize',6,'horiz','left','vert','middle','color',col(183,:),'fontweight','bold')
box off
axis tight
%ylabel 'change in GL flux (Gt/yr)'
ntitle('Response to Coastal Change','fontsize',8)
ylim([-2.5 19])
set(gca,'fontsize',7)


%export_fig observed_thinning_vs_calving_issm_response.png -r600 -p0.01
return
%%

figure
hold on
for k=1:181
   pl(k) = plot(0:100,glf_cf(k,:)-glf_cf(k,1),'-','color',col(k,:)); 
   if k<=181
      text(100,glf_cf(k,end)-glf_cf(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:))
   end
end
pl(183) = plot(0:100,glf_cf(183,:)-glf_cf(183,1),'-','color',col(183,:)); 
pl(183).LineWidth = 1; 
text(100,glf_cf(183,end)-glf_cf(183,1),'Antarctica','fontsize',6,'horiz','left','vert','middle','color',col(183,:),'fontweight','bold')
box off
axis tight
%ylabel 'change in GL flux (Gt/yr)'
ntitle('Response to Coastal Change','fontsize',8)
set(gca,'fontsize',7)


%%

xtmp = future_years;
ytmp = glf_tf-glf_tf(:,1);
figure
hold on
for k=1:181
   pl(k) = plot(xtmp,ytmp(k,:),'color',col(k,:)); 
   if k<=181
      text(xtmp(end),ytmp(k,end),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:))
   end
end
pl(183) = plot(xtmp,ytmp(183,:),'color',col(183,:)); 
pl(183).LineWidth = 2; 
text(xtmp(end),ytmp(183,end),'Antarctica','fontsize',6,'horiz','left','vert','middle','color',col(183,:),'fontweight','bold')
box off
axis tight
ylabel 'change in GL flux (Gt/yr)'
xlabel 'years into the future'

%%
M = load('calving_flux_timeseries.mat');

%figure 
%plot(M.year,M.M_calving(:,10))

figure
hs = scatter(M.M_calving(:,10),glf_cp(10,:),100,M.year,'filled');
xlabel 'ice shelf mass (Gt)' 
ylabel 'GL flux (Gt/yr)'
cb = colorbar; 

%%
k=k+1; 
hs.XData = M.M_calving(:,k); 
hs.YData = glf_cp(k,:); 
title(Dn.name{k})

%%


%%

flux_per_calve = nan(183,1); 
flux_per_calve_p = flux_per_calve; 
flux_per_calve_r = flux_per_calve; 
flux_per_calve_sig = false(size(flux_per_calve)); 

for k=1:183
   pv = polyfit(M.M_calving(:,k),glf_cp(k,:)',1); 
   flux_per_calve(k) = pv(1); 
   [flux_per_calve_r(k),flux_per_calve_p(k)] = corr(M.M_calving(:,k),glf_cp(k,:)'); 
   flux_per_calve_sig(k) = mann_kendall(M.M_calving(:,k),glf_cp(k,:)); 
end

%%

HM = load('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/hypothetical_iceshelf_mass.mat')

flux_per_calvef = nan(183,1); 
flux_per_calvef_sig = nan(183,1); 

ind = 1:51; 

for k=1:183
   pv = polyfit(HM.M_calving_future(ind,k),glf_cf(k,ind)',1); 
   flux_per_calvef(k) = 1000*pv(1); 
   flux_per_calvef_sig(k) = mann_kendall(HM.M_calving_future(ind,k),glf_cf(k,ind)'); 
end

figure
patchsc(Dn.x,Dn.y,flux_per_calvef(1:181),...
   'colormap',crameri('-vik'),...
   'caxis',[-1 1]*2,...
   'edgecolor','none')
axis tight off
bedmachine
cb=colorbarps; 
title 'calving sensitivity'

%%


HM = load('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/hypothetical_iceshelf_mass.mat')

flux_per_thinf = nan(183,1); 
flux_per_thinf_sig = nan(183,1); 

ind = 1:26; 

for k=1:183
   pv = polyfit(HM.M_thinning_future(ind,k),glf_tf(k,ind)',1); 
   flux_per_thinf(k) = 1000*pv(1); 
   flux_per_thinf_sig(k) = mann_kendall(HM.M_thinning_future(ind,k),glf_tf(k,ind)'); 
end

figure
patchsc(Dn.x,Dn.y,flux_per_thinf(1:181),...
   'colormap',crameri('-vik'),...
   'caxis',[-1 1]*2,...
   'edgecolor','none')
axis tight off
bedmachine
cb=colorbarps; 
title 'thinning sensitivity'

%%

ind = flux_per_thinf_sig&flux_per_calvef_sig;
flux_per_thinf(~ind) = nan; 

figure
scatter(flux_per_thinf(1:181),flux_per_calvef(1:181))

%axis([-1 1 -1 1]*150)

text(flux_per_thinf(1:181),flux_per_calvef(1:181),Dn.name,'horiz','center','vert','bot','fontsize',6)

%%

D = load('iceshelf_thickness_cube_1992-2018.mat'); 

%%

% Determine the mass of ice in each grid cell: 
[X,Y] = meshgrid(D.x,D.y); 
[Lat,~] = ps2ll(X,Y); 
F = 917*1e-12*(diff(D.x(1:2))./ psdistortion(Lat) ).^2; % F is a factor that converts vertical rates of ice change to Gt per grid cell.  
clear Lat X Y

F2 = cube2rect(F,D.always_ice & ~D.ground); 
H2 = double(cube2rect(D.H,D.always_ice & ~D.ground)); 
mask2 = cube2rect(permute(ncread(fn,'iceshelf_mask'),[2 1]),D.always_ice & ~D.ground); 

M_thinning = nan(26,183); 

for k = 1:181
   if any(mask2==k)
      M_thinning(:,k) = sum(H2(:,mask2==k).*F2(mask2==k),2); 
   end
end

M_thinning(:,182) = sum(H2(:,mask2==0).*F2(mask2==0),2); 
M_thinning(:,183) = sum(H2.*F2,2); 

clear D.H
%%

flux_per_thin = nan(183,1); 
flux_per_thin_p = flux_per_thin; 
flux_per_thin_r = flux_per_thin; 
flux_per_thin_sig = false(size(flux_per_thin)); 

for k=1:183
   pv = polyfit(M_thinning(:,k),glf_tp(k,:)',1); 
   flux_per_thin(k) = pv(1); 
   [flux_per_thin_r(k),flux_per_thin_p(k)] = corr(M_thinning(:,k),glf_tp(k,:)'); 
   flux_per_thin_sig(k) = mann_kendall(M_thinning(:,k),glf_tp(k,:)); 
end

%%

tmpt = flux_per_thin;
tmpc = flux_per_calve; 
tmpt(~flux_per_thin_sig) = 0; 
tmpc(~flux_per_calve_sig) = 0; 
tmpt(flux_per_thin_p>.001) = 0; 
tmpc(flux_per_calve_p>.001) = 0; 
tmpt(tmpt>0) = 0; 
tmpc(tmpc>0) = 0; 

tmp = tmpt; 
tmp(abs(tmpc)>abs(tmpt)) = -tmpc(abs(tmpc)>abs(tmpt)); 

figure
patchsc(Dn.x,Dn.y,1000*tmp(1:181),...
   'caxis',[-1 1]*2,...
   'colormap',crameri('-vik'),...
   'edgecolor','none')
axis tight off
bedmachine
cb=colorbarps; 


%%


figure('pos',[29.00        403.00        813.00        420.00])
subplot(1,2,1) 
patchsc(Dn.x,Dn.y,1000*flux_per_thin(1:181),...
   'caxis',[-1 1]*2,...
   'colormap',crameri('-vik'),...
   'edgecolor','none')
axis tight off
bedmachine
cb=colorbarps; 
xlabel(cb,'Mt/yr per Gt')
title 'sensitivity to thinning'

subplot(1,2,2) 
patchsc(Dn.x,Dn.y,1000*flux_per_calve(1:181),...
   'caxis',[-1 1]*2,...
   'colormap',crameri('-vik'),...
   'edgecolor','none')
axis tight off
bedmachine
cb=colorbarps; 
xlabel(cb,'Mt/yr per Gt')
title 'sensitivity to calving'

%%

k=10; 

figure
subplot(1,2,1)
hst = scatter(M_thinning(:,k)-M_thinning(1,k),glf_tp(k,:),100,D.year,'filled');
xlabel 'mass change due to thinning (Gt)' 
ylabel 'GL flux (Gt/yr)'
cb = colorbar('location','southoutside'); 
caxis([1992 2021])
hold on
pft = polyplot(M_thinning(:,k)-M_thinning(1,k),glf_tp(k,:),1,'r');
uistack(pft,'bottom')
ntt = ntitle([num2str(1000*flux_per_thin(k)),' Mt/yr per Gt'],'color','r');

subplot(1,2,2)
hsc = scatter(M.M_calving(:,k)-M.M_calving(1,k),glf_cp(k,:),100,M.year,'filled');
xlabel 'mass change due to calving (Gt)' 
ylabel 'GL flux (Gt/yr)'
cb = colorbar('location','southoutside'); 
caxis([1992 2021])
hold on
pfc = polyplot(M.M_calving(:,k)-M.M_calving(1,k),glf_cp(k,:),1,'r');
uistack(pfc,'bottom'); 
ntc = ntitle([num2str(1000*flux_per_calve(k)),' Mt/yr per Gt'],'color','r');

sgtitle(Dn.name{k})
k=0;
%%

k=k+1; 
delete(pft)
delete(pfc) 

subplot(1,2,1)
pft = polyplot(M_thinning(:,k)-M_thinning(1,k),glf_tp(k,:),1,'r');
uistack(pft,'bottom')
subplot(1,2,2) 
pfc = polyplot(M.M_calving(:,k)-M.M_calving(1,k),glf_cp(k,:),1,'r');
uistack(pfc,'bottom'); 

ntt.String = [num2str(1000*flux_per_thin(k)),' Mt/yr per Gt']; 
ntc.String = [num2str(1000*flux_per_calve(k)),' Mt/yr per Gt']; 

hst.XData = M_thinning(:,k)-M_thinning(1,k); 
hst.YData = glf_tp(k,:); 
hsc.XData = M.M_calving(:,k)-M.M_calving(1,k); 
hsc.YData = glf_cp(k,:); 
sgtitle(Dn.name{k})

%%

[~,lon] = ps2ll(Dn.x_center,Dn.y_center);
xtmp = 10*(glf_tp(:,end)-glf_tp(:,1))./(D.year(end)-D.year(1)); 
ytmp = 10*(glf_cp(:,end)-glf_cp(:,1))./(I.year(end)-I.year(1));


ind = hypot(xtmp(1:181),ytmp(1:181))>.08; 
figure
scatter(xtmp(1:181),ytmp(1:181),sqrt(Dn.area_km2/10),lon,'filled')
text(xtmp(ind),ytmp(ind),Dn.name(ind),'horiz','center','vert','bottom','fontsize',6)
hold on
axis tight
% pl=plot(xlim,xlim,'color',rgb('gray')); 
% uistack(pl,'bottom')
axis tight equal
box off
caxis([-1 1]*180)
cmocean phase
set(gca,'fontsize',7)
xlabel 'GL flux response to thickness change, 1992 to 2017 (Gt/decade)'
ylabel 'GL flux response to thickness change, 1997 to 2021 (Gt/decade)'

axis([-1 1 -1 1]*2.5)
ax_old = axis; 
axis(max(ax_old)*[-1 1 -1 1]); 

ax = axis; 
box on
%set(gca,'color',rgb('pastel pink')); 
cm = crameri('-vik',1001); % broc, cork, vik
set(gca,'color',cm(501-60,:)); 

hp = patch([ax(1) ax(2) 0 ax(1)],[ax(3) ax(3) 0 ax(3)],'b'); 
%hp.FaceColor = rgb('pastel blue'); 
hp.FaceColor = cm(501+60,:); 

hp.EdgeColor = 'none'; 
uistack(hp,'bottom') 

hp(2) = patch([ax(1) 0 ax(2) ax(1)],[ax(4) 0 ax(4) ax(4)],'b'); 
%hp(2).FaceColor = rgb('pastel blue'); 
hp(2).FaceColor = cm(501+60,:) ;

hp(2).EdgeColor = 'none'; 
uistack(hp(2),'bottom') 
axis(ax_old) 