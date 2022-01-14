
% Run this from issm_thickness_response_analysis.m.

%% Load data 

load('issm_gl_flux.mat'); % from issm_thickness_response_analysis 
Dn = load('iceshelves_2008_v2.mat');

if ~exist('P','var')
   for k = 1:181
      tmpx = Dn.x{k}; 
      tmpy = Dn.y{k}; 
      P(k) = polyshape(tmpx,tmpy); 
   end
end
%% Generate pdfs of response to thinning and calving. 

colt = hex2rgb('f9575a'); 
colc = hex2rgb('6f78f2'); 

Dn.name{182}='Other'; 
Dn.name{183}='Antarctica'; 

lw2 = 1; % composite linewidth
ms=3; % markersize

for kk=0:6
   figure('pos',[40 40 560 760])
   for k = 1:27
      if (kk*27+k)<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            plot(D.year,glf_tp(kk*27+k,:),'color',colt,'linewidth',lw2,'markersize',ms)
            hold on
            plot(I.year,glf_cp(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)

            box off
            axis tight
            %xlim([1997 2022])
            ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')

            ax = gca; 
            %ax.YAxis.Exponent = 0;
            %ytickformat('%.0f');
            set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])

      end
      if kk==6
         k=21; % just to set the legend on the last plot properly
      end
   end

   % Make a legend for the last axes: 
   lg = legend('response to thinning','response to calving','location','southwest');

   % Create a new axis just to get its position: 
   tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);

   % Move the legend to the empty axis position and delete the empty axes: 
   lg.Position=tmpax.Position;
   delete(tmpax)

   if kk==0
      sgtitle('Modeled grounding line flux (Gt yr^{-1})','fontsize',8) 
   else
      sgtitle('Modeled grounding line flux (Gt yr^{-1}), continued','fontsize',8) 
   end

   export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_thinning_and_calving_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')

end

%% Generate pdfs of response to hypothetical thinning

for kk=0:6
   figure('pos',[40 40 560 760])
   for k = 1:27
      if (kk*27+k)<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            plot(future_years,glf_tf(kk*27+k,:),'color',colt,'linewidth',lw2,'markersize',ms)
            hold on
            %plot(0:100,glf_cf(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)

            box off
            axis tight
            %xlim([1997 2022])
            ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')

            ax = gca; 
            %ax.YAxis.Exponent = 0;
            %ytickformat('%.0f');
            set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])

      end
      if kk==6
         k=21; % just to set the legend on the last plot properly
      end
   end

   % Make a legend for the last axes: 
   lg = legend('response to thinning','location','southwest');

   % Create a new axis just to get its position: 
   tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);

   % Move the legend to the empty axis position and delete the empty axes: 
   lg.Position=tmpax.Position;
   delete(tmpax)

   sgtitle('Modeled grounding line flux (Gt yr^{-1}) response to years of continued trends','fontsize',8) 

   export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_hypothetical_thinning_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')

end

%% Generate pdfs of response to hypothetical calving

for kk=0:6
   figure('pos',[40 40 560 760])
   for k = 1:27
      if (kk*27+k)<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            plot(0:100,glf_cf(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)
            hold on
            %
            box off
            axis tight
            %xlim([1997 2022])
            ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')

            ax = gca; 
            %ax.YAxis.Exponent = 0;
            %ytickformat('%.0f');
            set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])

      end
      if kk==6
         k=21; % just to set the legend on the last plot properly
      end
   end

   % Make a legend for the last axes: 
   lg = legend('response to calving','location','southwest');

   % Create a new axis just to get its position: 
   tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);

   % Move the legend to the empty axis position and delete the empty axes: 
   lg.Position=tmpax.Position;
   delete(tmpax)

   sgtitle('Modeled grounding line flux (Gt yr^{-1}) response to ice shelf loss, by percent area','fontsize',8) 

   export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_hypothetical_calving_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')

end

%% Total ice shelf collapse scenario 

% Change in GL flux (percent): 
dg = 100*(glf_cf(:,end)-glf_cf(:,1))./glf_cf(:,1);

Asc = sqrt(Dn.area_km2); 
Asc = 90*Asc/max(Asc); 

[~,lon] = ps2ll(Dn.x_center,Dn.y_center);

figure('pos',[200 200 308 231])
scatter(glf_cf(1:181,1),dg(1:181),Asc,lon,'filled')
axis tight
box off
caxis([-180 180]) 
cmocean phase 
ind = glf_cf(:,end)>25; 
ind(dg>200) = true; 
ind([19 33 49 74 81 87 113 146 147 182:183]) = false; % remove labels that clutter 
text(glf_cf(ind,1),dg(ind),Dn.name(ind),'horiz','center','vert','bot','fontsize',6,...
   'fontangle','italic','color',.4*[1 1 1])
axis tight
title('Total ice shelf collapse scenario','fontsize',8)
xlabel('Modeled present-day GL flux (Gt/yr)','fontsize',7)
ylabel('Percent change in GL flux','fontsize',7)
set(gca,'fontsize',7)

[X,Y] = meshgrid(0:10:250,0:10:450);

%G = (1+Y/100).*X; 
G = (Y/100).*X; 
hold on

[C,hC] = contour(X,Y,G,0:50:1500,'k'); 
hC.Color = .8*[1 1 1];
hC.LineWidth = 0.3; 

uistack(hC,'bottom')
axis([0 250 0 450])

n = 50:50:1100; 
ny = 100*n./(250*ones(size(n))); 
str = num2str(n','%-4.f'); 
txt = text(250*ones(size(n)),ny,str,'vert','middle','fontsize',6,'color',.8*[1 1 1]); 
col = mat2rgb(lon,cmocean('phase'),[-1 1]*180); 

gp = plotboxpos(gca);
axes('position',[gp(1)+(.75)*gp(3) gp(2) .25*gp(3) .25*gp(4)])
 
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
set(gca,'xcolor','none','ycolor','none','pos',[gp(1)+(.75)*gp(3) gp(2) gp2(3) gp2(4)])
set(gca,'fontsize',7)

% export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/catastrophic_calving_response.jpg','-pdf','-r600','-painters')
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
