% This script loads and plots data from the extruded_antarctica*.h5 file
% that has already been created by flow_dem_extend.m. This is just to
% verify that everything went okay with the extrusion process. 
% 
% Chad A. Greene, NASA JPL, Sept 2021. 

%% Inspect data: 

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-09-01.h5'; 
h5disp(fn)

%% Load data: 

x = h5read(fn,'/x'); 
y = h5read(fn,'/y'); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 
v_source = permute(h5read(fn,'/v_source'),[2 1]); 

id = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
H_source = permute(h5read(fn,'/thickness_source'),[2 1]); 

%% Ice Shelf mask 

figure
imagescn(x,y,id)
bedmachine
axis off 
cm = rand(182,3); 
cm(1,:) = [1 1 1];
colormap(cm)
title('iceshelf_mask','interpreter','none','fontname','courier')

% export_fig extruded_iceshelf_mask.png -r500 -p0.01 

%% Velocity Source 

figure
imagescn(x,y,v_source) 
bedmachine
axis off
ax(1) = gca; 
% cm = hex2rgb({'ffffff';'ad993c';'7066bc';'56ae6c';'b54f90';'b94e45'});
% colormap(cm)

% Categorical colormap, sorted by lightness to reflect prioritization 
cm = brewermap(5,'set1'); 
lab = rgb2lab(cm); 
[~,idx] = sort(lab(:,1),'desc'); 
cm = cm(idx,:); 
colormap(cm); 

caxis([0.5 5.5])
cb = colorbar; 
set(cb,'ytick',1:5,'yticklabel',{'ITSLIVE';'MEaSUREs v2';'blend';'interpolated';'extrapolated'},'ydir','reverse')
title('iceshelf_mask','interpreter','none','fontname','courier')

% export_fig extruded_v_source.png -r500 -p0.01 

%% Velocity 

figure
imagescn(x,y,hypot(vx,vy))
cmocean thermal
hgl = bedmachine('gl','color','k');
hold on
axis off
set(gca,'colorscale','log')
caxis([1.5 4000])

% Optional quiver plot: 
if false
   hold on
   q = quiversc(x,y,vx,vy,'density',200);
   q.Color = rgb('dark red');
   q.AutoScaleFactor = 2;
end
title('velocity','interpreter','none','fontname','courier')
cb = colorbar; 
ylabel(cb,'velocity m/yr') 

linkaxes(ax,'xy') 

% export_fig extruded_v.png -r500 -p0.01 

%% Thickness Source 

figure
imagesc(x,y,H_source)
axis xy off
bedmachine 

% Categorical colormap, sorted by lightness to reflect prioritization 
cm = brewermap(7,'accent'); 
lab = rgb2lab(cm); 
[~,idx] = sort(lab(:,1),'desc'); 
cm = cm(idx,:); 
colormap(cm); 

title('thickness_source','interpreter','none','fontname','courier')
caxis([0.5 7.5])
cb = colorbar; 
set(cb,'ytick',1:7,'yticklabel',{'BedMachine';'REMA-GEMB';'Bedmap2-GEMB';'Bamber-GEMB';'RAMP2-GEMB';'interpolated';'extrapolated'},'ydir','reverse')

ax(1) = gca; 

% export_fig extruded_thickness_source.png -r500 -p0.01 

%% Thickness 

figure
imagescn(x,y,H)
axis xy off
bedmachine 
ax(2) = gca; 
caxis([0 500]) 
title('thickness','interpreter','none','fontname','courier')

linkaxes(ax,'xy') 

% export_fig extruded_thickness.png -r500 -p0.01 

