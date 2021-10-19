% This script loads and plots data from the extruded_antarctica*.h5 file
% that has already been created by flow_dem_extend.m. This is just to
% verify that everything went okay with the extrusion process. 
% 
% Chad A. Greene, NASA JPL, Sept 2021. 

%% Inspect data: 

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
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
bedmachine('coast','color',rgb('gray'),'linewidth',0.25)
bedmachine('gl','color',rgb('gray'),'linewidth',0.25)
axis off 
cm = rand(182,3); 
cm(1,:) = [1 1 1];
colormap(cm)
ntitle('iceshelf_mask','interpreter','none','fontname','courier','fontweight','bold','fontsize',14,'location','center')
set(gca,'pos',[0 0 1 1])

% export_fig extruded_iceshelf_mask.png -r500 -painters -nocrop

%% Velocity Source 

v_source(v_source==0) = 5; 
figure
h = imagescn(x,y,v_source) ; 
bedmachine('coast','color',rgb('gray'),'linewidth',0.25)
bedmachine('gl','color',rgb('gray'),'linewidth',0.25)
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
set(gca,'pos',[0 0 1 1])

cb = colorbarps; 
set(cb,'xtick',1:5,'xticklabel',{'ITSLIVE';'MEaSUREs v2';'blend';'interpolated';'extrapolated'})
ntitle('v_source','interpreter','none','fontname','courier','fontweight','bold','fontsize',14,'location','center')

% export_fig extruded_v_source.png -r500 -painters -nocrop 

%% Velocity 

figure
imagescn(x,y,hypot(vx,vy))
cmocean thermal
hgl = bedmachine('gl','color','k','linewidth',0.25);
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
set(gca,'pos',[0 0 1 1])
ntitle('velocity','interpreter','none','fontname','courier','fontweight','bold','fontsize',14,'color','w','location','center')
cb = colorbarps; 
xlabel(cb,'velocity m/yr') 
set(cb,'xtick',[10 100 1000],'color','w')

linkaxes(ax,'xy') 

% export_fig extruded_velocity.png -r500 -painters -nocrop 

%% Thickness Source 

figure
imagesc(x,y,H_source)
axis xy off
bedmachine('coast','color',rgb('gray'),'linewidth',0.25)
bedmachine('gl','color',rgb('gray'),'linewidth',0.25)

% Categorical colormap, sorted by lightness to reflect prioritization 
cm = brewermap(7,'accent'); 
lab = rgb2lab(cm); 
[~,idx] = sort(lab(:,1),'desc'); 
cm = cm(idx,:); 
colormap(cm); 

set(gca,'pos',[0 0 1 1])
ntitle('thickness_source','interpreter','none','fontname','courier','fontweight','bold','fontsize',14,'location','center')
caxis([0.5 7.5])
cb = colorbarps; 
set(cb,'xtick',1:7,'xticklabel',{'BedMachine';'REMA-GEMB';'Bedmap2-GEMB';'Bamber-GEMB';'RAMP2-GEMB';'interpolated';'extrapolated'})

ax(1) = gca; 

% export_fig extruded_thickness_source.png -r500 -painters -nocrop 

%% Thickness 

figure
imagescn(x,y,H)
axis xy off
bedmachine('coast','color',rgb('gray'),'linewidth',0.25)
bedmachine('gl','color',rgb('gray'),'linewidth',0.25)
ax(2) = gca; 
caxis([0 500]) 
set(gca,'pos',[0 0 1 1])
ntitle('thickness','interpreter','none','fontname','courier','fontweight','bold','fontsize',14,'color','w','location','center')
cmocean dense
linkaxes(ax,'xy') 
cb=colorbarps; 
xlabel(cb,'thickness (m)')
% export_fig extruded_thickness.png -r500 -painters -nocrop 

