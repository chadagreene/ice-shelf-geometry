% This script converts the icemask_composite.mat ice cube into a timeseries
% of calving flux. 
% 
% Chad Greene October 2021. 
% NASA Jet Propulsion Laboratory. 


%% Load data: 

load icemask_composite.mat

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

[X,Y] = meshgrid(x,y); 
rho_ice = 917;
D = load('iceshelves_2008_v2.mat');
tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]); 

%%

always_ice = imfill(all(ice,3),8,'holes'); 
ever_ice = imfill(any(ice,3),8,'holes'); 
action_ice = ever_ice & ~always_ice; % the grid cells where anything ever changes

% Determine the mass of ice in each grid cell: 
[Lat,~] = ps2ll(X,Y); 
A = (diff(x(1:2)) ./ psdistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  
IceMass = rho_ice.*H.*A.*1e-12; % ice mass in each grid cell (Gt)
clear Lat

% Position of grid cells after dt: 
dt = 0.001; % one thousandth of a year. This value doesn't affect the final results because they're multiplied by dt/dt. Just ensure that dt is small enough that the fastest-moving ice doesn't move more than half a pixel or so over the course of dt.    
X1 = X + vx.*dt; % dt is very small, so vx.*dt will probably only be a few meters.
Y1 = Y + vy.*dt; 

%% Steady-state mass rate: 

% A grid of steady-state mass rate: 
Mdot_ss = (interp2(x,y,double(~always_ice),X1,Y1) - double(~always_ice)) .* IceMass./dt; 

% Convert the grid to 2D for quick summing: 
Mdot_ss2 = cube2rect(Mdot_ss,ever_ice & isfinite(Mdot_ss)); 
mask2 = cube2rect(mask,ever_ice & isfinite(Mdot_ss)); 

CR_ss = nan(183,1); 
for k = 1:181
   CR_ss(k) = sum(Mdot_ss2(:,mask2==k)); 
end
CR_ss(182) = sum(Mdot_ss2(:,mask2==0)); 
CR_ss(183) = sum(Mdot_ss2);
CR_ss(183)

figure
bar((1:182)-.5,CR_ss(1:182))
set(gca,'xtick',.5:181.5,'xticklabel',D.name)
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

%%

col = mat2rgb(CR_ss(1:181),cmocean('amp')); 
figure
hold on
p = plot(P,'edgecolor','none','facealpha',0.8); 
for k = 1:181
   p(k).FaceColor = col(k,:); 
end

bedmachine('gl','color',rgb('gray'))
axis tight off 

%% Change in mass

ice2 = cube2rect(ice,action_ice); 
IceMass2 = cube2rect(IceMass,action_ice); 
mask2 = cube2rect(mask,action_ice); 

% Using a negative sign here so positive is calving (to match steady state
% massrate) 
M2 = double(ice2).*IceMass2;

M_anom = nan(24,183); 
for k = 1:181
   M_anom(:,k) = sum(M2(:,mask2==k),2); 
end
M_anom(:,182) = sum(M2(:,mask2==0),2); 
M_anom(:,183) = sum(M2,2); 

dM = diff(M_anom); 


figure
hold on
for k = 1:23
   plot([year(k) year(k+1)],CR_ss(end) + dM(k,end)./diff([year(k) year(k+1)]).*[1 1],'linewidth',2,'color',rgb('dark blue'))
end
box off
axis tight
yline(CR_ss(end))
ylabel 'calving rate (Gt/yr)' 

%% 

Mtmp = CR_ss - ((M_anom(end,:) - M_anom(1,:))/(year(end)-year(1)))';

ind = CR_ss>1 | Mtmp>1; 
ind(182:183) = false; 
f=find(ind); 

figure
plot(CR_ss(1:181),Mtmp(1:181)','.')
for k = 1:length(f)
   txt(k) = text(CR_ss(f(k)),Mtmp(f(k)),D.name{f(k)},...
   'vert','bot','fontsize',ceil(abs((Mtmp(f(k))-CR_ss(f(k))))/5),'horiz','center');
end
axis image 
hold on
h = plot(xlim,xlim,'color',rgb('gray')); 
uistack(h,'bottom')
box off
xlabel 'steady-state calving rate (Gt/yr)'
ylabel 'observed mean calving rate 1997-2021 (Gt/yr)' 

% export_fig '/Users/cgreene/Documents/MATLAB/DEM_generation/calving_steadystate_vs_observed.png' -r500 -p0.01 -painters

%%

Mtmp = (M_anom(end,:) - M_anom(1,:))' - CR_ss.*(year(end)-year(1));


col = mat2rgb(100*(Mtmp(1:181)-CR_ss(1:181))./CR_ss(1:181),cmocean('bal'),[-1 1]*50); 
figure
hold on
p = plot(P,'edgecolor','none','facealpha',1); 
for k = 1:181
   p(k).FaceColor = col(k,:); 
end

bedmachine('gl','color',rgb('gray'))
axis tight off 

cb = colorbarps;
xlabel(cb,'Calving mass balance (Gt/yr)')
caxis([-1 1]*50)
cmocean bal

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

