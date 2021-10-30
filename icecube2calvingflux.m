% This script converts the icemask_composite.mat ice cube into a timeseries
% of calving flux. 
% 
% Creates calving_flux_timeseries.mat
% 
% Chad Greene October 2021. 
% NASA Jet Propulsion Laboratory. 

%% Load data: 

load icemask_composite.mat
load('issm_calving_melt_setup.mat','ground')

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

[X,Y] = meshgrid(x,y); 
rho_ice = 917;
D = load('iceshelves_2008_v2.mat');
tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]) & ~ground; 

%%

always_ice = imfill(all(ice,3),8,'holes'); 
ever_ice = imfill(any(ice,3),8,'holes'); 
action_ice = ever_ice & ~always_ice; % the grid cells where anything ever changes

% Determine the mass of ice in each grid cell: 
[Lat,~] = ps2ll(X,Y); 
A = (diff(x(1:2)) ./ psdistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  
IceMass = rho_ice.*H.*A.*1e-12; % ice mass in each grid cell (Gt)
clear Lat

%% Position of grid cells after one year
% It is possible to calculate the displacements after a very tiny dt value and then divide by 
% that small dt, however, when I did this in a test against Alex's flux gate method, I noticed 
% some numerical noise along all the near-zero-velocity pixels of ice affected the continent-wide 
% GL flux when very small dt values are used. Therefore, I think it's easier and physically justified 
% to use pixel locations after 1 year of displacement. But to account for acceleration and curves, 
% use several sub-year steps to get to 1 year. 

X1 = X; 
Y1 = Y; 

for k = 1:10
   dx = interp2(x,y,vx,X1,Y1); 
   dy = interp2(x,y,vy,X1,Y1); 
   
   X1 = X1+dx/10; 
   Y1 = Y1+dy/10; 
   k
end

%% Steady-state mass rate: 

% A grid of steady-state mass rate at the calving front: 
Mdot_cf = (double(always_ice) - interp2(x,y,double(always_ice),X1,Y1)) .* IceMass; 
Mdot_cf(isnan(Mdot_cf)) = 0; % bc we count it as zero anyway, and now we don't have to deal with nans. 

% A grid of steady-state mass rate at the calving front: 
Mdot_gl = (interp2(x,y,double(tidal),X1,Y1) - double(tidal)) .* IceMass; 
Mdot_gl(isnan(Mdot_gl)) = 0; 

% Convert the grid to 2D for quick summing: 
Mdot_cf2 = cube2rect(Mdot_cf,ever_ice); 
Mdot_gl2 = cube2rect(Mdot_gl,ever_ice); 
mask2 = cube2rect(mask,ever_ice); 

CF_ss = nan(183,1); 
GL_ss = nan(183,1); 
for k = 1:181
   CF_ss(k) = sum(Mdot_cf2(:,mask2==k)); 
   GL_ss(k) = sum(Mdot_gl2(:,mask2==k)); 
end
CF_ss(182) = sum(Mdot_cf2(:,mask2==0)); 
CF_ss(183) = sum(Mdot_cf2);
CF_ss(183)

GL_ss(182) = sum(Mdot_gl2(:,mask2==0)); 
GL_ss(183) = sum(Mdot_gl2);
GL_ss(183)

figure
bar((1:182)-.5,CF_ss(1:182))
hold on
bar((1:182)-.5,-GL_ss(1:182))
set(gca,'xtick',.5:181.5,'xticklabel',D.name)
box off
axis tight

figure
plot(GL_ss(1:181),CF_ss(1:181),'.') 
hold on
text(GL_ss(1:181),CF_ss(1:181),D.name,'horiz','center','vert','bot','fontsize',6,'clipping','on')
box off
axis tight

%%

TotalAreaChange = (sum(A(ice(:,:,end)),'all') - sum(A(ice(:,:,1)),'all'))/(1000^2)
TotalMassChange = sum(IceMass(ice(:,:,end)),'all') - sum(IceMass(ice(:,:,1)),'all')

%% Map of steady state calving rates

warning off
for k = 1:181
   P(k) = polyshape(D.x{k},D.y{k}); 
end

col = mat2rgb(CF_ss(1:181),cmocean('amp')); 
figure
hold on
p = plot(P,'edgecolor','none','facealpha',0.8); 
for k = 1:181
   p(k).FaceColor = col(k,:); 
end

bedmachine('gl','color',rgb('gray'))
axis tight off 

%% Area and Mass time series

ice2 = cube2rect(ice,tidal & ever_ice); 
IceMass2 = cube2rect(IceMass,tidal & ever_ice); 
mask2 = cube2rect(mask,tidal & ever_ice); 

M2 = double(ice2).*IceMass2;
A2 = double(ice2).*cube2rect(A,tidal & ever_ice);

M_calving = nan(24,183); 
A_calving = nan(24,183); 
for k = 1:181
   M_calving(:,k) = sum(M2(:,mask2==k),2); 
   A_calving(:,k) = sum(A2(:,mask2==k),2); 
end
M_calving(:,182) = sum(M2(:,mask2==0),2); 
M_calving(:,183) = sum(M2,2); 
A_calving(:,182) = sum(A2(:,mask2==0),2); 
A_calving(:,183) = sum(A2,2); 

dM = diff(M_calving); 

figure
hold on
for k = 1:23
   plot([year(k) year(k+1)],CF_ss(end) - dM(k,end)./diff([year(k) year(k+1)]).*[1 1],'linewidth',2,'color',rgb('dark blue'))
end
box off
axis tight
yline(CF_ss(end))
ylabel 'calving rate (Gt/yr)' 

readme = 'ice area time series and corresponding mass time series.  Created by icecube2calvingflux.m';   
% save('calving_flux_timeseries.mat','year','readme','CF_ss','GL_ss','A_calving','M_calving')

%% 
% 
% Mtmp = CF_ss - ((M_calving(end,:) - M_calving(1,:))/(year(end)-year(1)))';
% 
% ind = CF_ss>1 | Mtmp>1; 
% ind(182:183) = false; 
% f=find(ind); 
% 
% figure
% plot(CF_ss(1:181),Mtmp(1:181)','.')
% for k = 1:length(f)
%    txt(k) = text(CF_ss(f(k)),Mtmp(f(k)),D.name{f(k)},...
%    'vert','bot','fontsize',ceil(abs((Mtmp(f(k))-CF_ss(f(k))))/5),'horiz','center');
% end
% axis image 
% hold on
% h = plot(xlim,xlim,'color',rgb('gray')); 
% uistack(h,'bottom')
% box off
% xlabel 'steady-state calving rate (Gt/yr)'
% ylabel 'observed mean calving rate 1997-2021 (Gt/yr)' 
% 
% % export_fig '/Users/cgreene/Documents/MATLAB/DEM_generation/calving_steadystate_vs_observed.png' -r500 -p0.01 -painters
% 
% %%
% 
% Mtmp = (M_calving(end,:) - M_calving(1,:))' - CF_ss.*(year(end)-year(1));
% 
% 
% col = mat2rgb(100*(Mtmp(1:181)-CF_ss(1:181))./CF_ss(1:181),cmocean('bal'),[-1 1]*50); 
% figure
% hold on
% p = plot(P,'edgecolor','none','facealpha',1); 
% for k = 1:181
%    p(k).FaceColor = col(k,:); 
% end
% 
% bedmachine('gl','color',rgb('gray'))
% axis tight off 
% 
% cb = colorbarps;
% xlabel(cb,'Calving mass balance (Gt/yr)')
% caxis([-1 1]*50)
% cmocean bal
% 
% %%


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

