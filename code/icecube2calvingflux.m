% This script converts the icemask_composite.mat ice cube into a timeseries
% of calving flux. 
% 
% Creates calving_flux_timeseries.mat
% 
% This script takes about 15 minutes to run on my laptop. 
% 
% Chad Greene October 2021. 
% NASA Jet Propulsion Laboratory. 

%% Load data: 

load icemask_composite.mat
load('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','ground')

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

[X,Y] = meshgrid(x,y); 
rho_ice = 917;
D = load('iceshelves_2008_v2.mat');
tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]) & ~ground; 

%% Error fields  

v_source = permute(h5read(fn,'/v_source'),[2 1]); 

vx_err = itslive_data('vx_err'); 
vy_err = itslive_data('vy_err'); 

% Use Measures v2 errors wherever it was the velocity source:
filename = 'antarctica_ice_velocity_450m_v2.nc'; 
xm = -2800000:450:2800000; 
ym = (-2800000:450:2800000)'; 
vx_err(v_source==2) = interp2(xm,ym,double(rot90(ncread(filename,'ERRX'))),X(v_source==2),Y(v_source==2)); 
vy_err(v_source==2) = interp2(xm,ym,double(rot90(ncread(filename,'ERRY'))),X(v_source==2),Y(v_source==2)); 

% Extend, but downsample first to make the inpainting tractable 
skip = 4; 
vx_err_r = vx_err(1:skip:end,1:skip:end); 
vy_err_r = vy_err(1:skip:end,1:skip:end); 
vx_err_r = regionfill(vx_err_r,isnan(vx_err_r)); 
vy_err_r = regionfill(vy_err_r,isnan(vy_err_r)); 
vx_err(isnan(vx_err)) = interp2(x(1:skip:end),y(1:skip:end),vx_err_r,X(isnan(vx_err)),Y(isnan(vy_err)));
vy_err(isnan(vy_err)) = interp2(x(1:skip:end),y(1:skip:end),vy_err_r,X(isnan(vy_err)),Y(isnan(vy_err)));

clear v_source filename xm ym vx_err_r vy_err_r skip
%%

always_ice = imfill(all(ice,3),8,'holes'); 
ever_ice = imfill(any(ice,3),8,'holes'); 
action_ice = ever_ice & ~always_ice; % the grid cells where anything ever changes

% Determine the mass of ice in each grid cell: 
[Lat,~] = ps2ll(X,Y); 
A = (diff(x(1:2)) ./ psdistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  
IceMass = rho_ice.*H.*A.*1e-12; % ice mass in each grid cell (Gt)
IceMass_plus = rho_ice.*(H+10).*A.*1e-12; % Thickness is primarily a hydrostatic inversion applied to REMA. REMA errors are not explicitly given for the mosaic, but the Cryosphere mansucript describing REMA says surface elevation errors are "less than 1 m" and the pdf description of the accuracy of the strips says errors average around 0.6 m. Hydrostatic inversion applied to 0.6 would be 5.6 m thickness uncertainty, however thickness uncertainty is also affected by FAC uncertainty. Now we're into uncertainties of uncertainties and how they compound together and we could put a great deal of effort into getting an exact measure of uncertainty, but of course even that would all be made up. Thus, let's just say there's 10 m thickness uncertainty everywhere and leave it at that.  
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

% Do it again with velocity error: 
X1_plus = X; 
Y1_plus = Y; 
for k = 1:10
   dx = interp2(x,y,(abs(vx)+vx_err).*sign(vx),X1_plus,Y1_plus); 
   dy = interp2(x,y,(abs(vy)+vy_err).*sign(vy),X1_plus,Y1_plus); 
   
   X1_plus = X1_plus+dx/10; 
   Y1_plus = Y1_plus+dy/10; 
   k
end

clear dx dy 
%% Steady-state mass rate: 

% A grid of steady-state mass rate at the calving front: 
Mdot_cf = (double(always_ice) - interp2(x,y,double(always_ice),X1,Y1)) .* IceMass; 
Mdot_cf(isnan(Mdot_cf)) = 0; % bc we count it as zero anyway, and now we don't have to deal with nans. 

% Error estimates: Note these are the multiplied effects of velocity and thickness error. Technically, vel and thicness errors are independent and their effects should be root-sum-squared, so these numbers may err on the side of being higher than formal errors.  
Mdot_cf_plus = (double(always_ice) - interp2(x,y,double(always_ice),X1_plus,Y1_plus)) .* IceMass_plus; 
Mdot_cf_plus(isnan(Mdot_cf_plus)) = 0; 

% A grid of steady-state mass rate at the calving front: 
Mdot_gl = (interp2(x,y,double(tidal),X1,Y1) - double(tidal)) .* IceMass; 
Mdot_gl(isnan(Mdot_gl)) = 0; 

Mdot_gl_plus = (interp2(x,y,double(tidal),X1_plus,Y1_plus) - double(tidal)) .* IceMass_plus; 
Mdot_gl_plus(isnan(Mdot_gl_plus)) = 0; 

% Convert the grid to 2D for quick summing: 
Mdot_cf2 = cube2rect(Mdot_cf,ever_ice); 
Mdot_cf2_plus = cube2rect(Mdot_cf_plus,ever_ice); 
Mdot_gl2 = cube2rect(Mdot_gl,ever_ice); 
Mdot_gl2_plus = cube2rect(Mdot_gl_plus,ever_ice); 
mask2 = cube2rect(mask,ever_ice); 

% Calculte the calving flux of each ice shelf: 
CF_ss = nan(183,1); 
CF_ss_plus = nan(183,1); 
GL_ss = nan(183,1); 
GL_ss_plus = nan(183,1); 
for k = 1:181
   CF_ss(k) = sum(Mdot_cf2(:,mask2==k)); 
   CF_ss_plus(k) = sum(Mdot_cf2_plus(:,mask2==k)); 
   GL_ss(k) = sum(Mdot_gl2(:,mask2==k)); 
   GL_ss_plus(k) = sum(Mdot_gl2_plus(:,mask2==k)); 
end
CF_ss(182) = sum(Mdot_cf2(:,mask2==0)); 
CF_ss(183) = sum(Mdot_cf2);
CF_ss_plus(182) = sum(Mdot_cf2_plus(:,mask2==0)); 
CF_ss_plus(183) = sum(Mdot_cf2_plus);

GL_ss(182) = sum(Mdot_gl2(:,mask2==0)); 
GL_ss(183) = sum(Mdot_gl2);
GL_ss_plus(182) = sum(Mdot_gl2_plus(:,mask2==0)); 
GL_ss_plus(183) = sum(Mdot_gl2_plus);

% Calculate error estimates as "plus" estimates minus 
CF_ss_err = CF_ss_plus - CF_ss; 
GL_ss_err = GL_ss_plus - GL_ss; 


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

if false 
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
end

%% Area and Mass time series

ice2 = cube2rect(ice,tidal & ever_ice); 
IceMass2 = cube2rect(IceMass,tidal & ever_ice); 
IceMass2_plus = cube2rect(IceMass_plus,tidal & ever_ice); 
mask2 = cube2rect(mask,tidal & ever_ice); 

M2 = double(ice2).*IceMass2;
M2_plus = double(ice2).*IceMass2_plus;
A2 = double(ice2).*cube2rect(A,tidal & ever_ice);

% We're focused on the floating pixels that change, but to account for the whole continent we'll need to add the grounded ice pixels too:  
A_neglected = sum(A(ever_ice & ~tidal),'all'); 
M_neglected = sum(IceMass(ever_ice & ~tidal),'all'); 

M_calving = nan(24,183); 
M_calving_plus = nan(24,183); 
A_calving = nan(24,183); 
for k = 1:181
   M_calving(:,k) = sum(M2(:,mask2==k),2); 
   M_calving_plus(:,k) = sum(M2_plus(:,mask2==k),2); 
   A_calving(:,k) = sum(A2(:,mask2==k),2); 
end
M_calving(:,182) = sum(M2(:,mask2==0),2)+M_neglected; 
M_calving(:,183) = sum(M2,2)+M_neglected; 
M_calving_plus(:,182) = sum(M2_plus(:,mask2==0),2)+M_neglected; 
M_calving_plus(:,183) = sum(M2_plus,2)+M_neglected; 
A_calving(:,182) = sum(A2(:,mask2==0),2)+A_neglected; 
A_calving(:,183) = sum(A2,2)+A_neglected; 

M_calving_err = M_calving_plus - M_calving; 

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
% save('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/calving_flux_timeseries.mat','year','readme','CF_ss','GL_ss','A_calving','M_calving','CF_ss_err','GL_ss_err','M_calving_err')

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

