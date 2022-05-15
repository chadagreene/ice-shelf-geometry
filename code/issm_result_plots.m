
% Run this from issm_thickness_response_analysis.m.

%% Load data 

M = load('calving_flux_timeseries.mat');
HM = load('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/hypothetical_iceshelf_mass.mat');
D = load('/Users/cgreene/Documents/MATLAB/DEM_generation/iceshelf_thickness_cube_1992-2018.mat','x','y','year');
I = load('icemask_composite.mat','year');
load('issm_gl_flux_strict.mat'); % from issm_thickness_response_analysis 
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

figure('pos',[40 40 560 760])

for kk=0:6
   clf
   for k = 1:27
      if (kk*27+k)<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            plot(D.year,glf_tp(kk*27+k,:),'color',colt,'linewidth',lw2,'markersize',ms)
            hold on
            plot(I.year,glf_cp(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)

            box off
            axis tight
            ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')

            ax = gca; 
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

   sgtitle({'Modeled grounding line flux (Gt yr^{-1})';'Instantaneous response to observed forcing'},'fontsize',8) 

%   export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_thinning_and_calving_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')
   filename = '/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_instantaneous_responses.pdf';
   if k==0
      assert(~exist(filename,'file'),[filename,' already exists'])
   end
   export_fig(filename,'-r600','-painters','-nocrop','-append')

end

%% Hypothetical ice shelf loss 

for kk=0:6
   clf
   for k = 1:27
      if (kk*27+k)<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            plot(0:100,glf_tf2(kk*27+k,:),'-','color',colt,'linewidth',lw2,'markersize',ms)
            hold on
            plot(0:100,glf_cf(kk*27+k,:),'-','color',colc,'linewidth',lw2,'markersize',ms)

            box off
            axis tight
            ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')

            ax = gca; 
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

   sgtitle({'Modeled grounding line flux (Gt yr^{-1})';'Instantaneous response to hypothetical percent mass loss'},'fontsize',8) 

   export_fig(filename,'-r600','-painters','-nocrop','-append')

end

return

%%

% Change in GL flux (percent): 
dg = 100*(glf_cf(:,end)-glf_cf(:,1))./glf_cf(:,1);

Asc = sqrt(Dn.area_km2); 
Asc = 90*Asc/max(Asc); 

[~,lon] = ps2ll(Dn.x_center,Dn.y_center);

% col = hex2rgb({'#000000','#FFFF00','#1CE6FF','#FF34FF','#FF4A46','#008941','#006FA6','#A30059','#FFDBE5','#7A4900','#0000A6','#63FFAC','#B79762','#004D43','#8FB0FF','#997D87','#5A0007','#809693','#FEFFE6','#1B4400','#4FC601','#3B5DFF','#4A3B53','#FF2F80','#61615A','#BA0900','#6B7900','#00C2A0','#FFAA92','#FF90C9','#B903AA','#D16100','#DDEFFF','#000035','#7B4F4B','#A1C299','#300018','#0AA6D8','#013349','#00846F','#372101','#FFB500','#C2FFED','#A079BF','#CC0744','#C0B9B2','#C2FF99','#001E09','#00489C','#6F0062','#0CBD66','#EEC3FF','#456D75','#B77B68','#7A87A1','#788D66','#885578','#FAD09F','#FF8A9A','#D157A0','#BEC459','#456648','#0086ED','#886F4C','#34362D','#B4A8BD','#00A6AA','#452C2C','#636375','#A3C8C9','#FF913F','#938A81','#575329','#00FECF','#B05B6F','#8CD0FF','#3B9700','#04F757','#C8A1A1','#1E6E00','#7900D7','#A77500','#6367A9','#A05837','#6B002C','#772600','#D790FF','#9B9700','#549E79','#FFF69F','#201625','#72418F','#BC23FF','#99ADC0','#3A2465','#922329','#5B4534','#FDE8DC','#404E55','#0089A3','#CB7E98','#A4E804','#324E72','#6A3A4C','#83AB58','#001C1E','#D1F7CE','#004B28','#C8D0F6','#A3A489','#806C66','#222800','#BF5650','#E83000','#66796D','#DA007C','#FF1A59','#8ADBB4','#1E0200','#5B4E51','#C895C5','#320033','#FF6832','#66E1D3','#CFCDAC','#D0AC94','#7ED379','#012C58','#7A7BFF','#D68E01','#353339','#78AFA1','#FEB2C6','#75797C','#837393','#943A4D','#B5F4FF','#D2DCD5','#9556BD','#6A714A','#001325','#02525F','#0AA3F7','#E98176','#DBD5DD','#5EBCD1','#3D4F44','#7E6405','#02684E','#962B75','#8D8546','#9695C5','#E773CE','#D86A78','#3E89BE','#CA834E','#518A87','#5B113C','#55813B','#E704C4','#00005F','#A97399','#4B8160','#59738A','#FF5DA7','#F7C9BF','#643127','#513A01','#6B94AA','#51A058','#A45B02','#1D1702','#E20027','#E7AB63','#4C6001','#9C6966','#64547B','#97979E','#006A66','#391406','#F4D749','#0045D2','#006C31','#DDB6D0','#7C6571','#9FB2A4','#00D891','#15A08A','#BC65E9','#FFFFFE','#C6DC99','#203B3C','#671190','#6B3A64','#F5E1FF','#FFA0F2','#CCAA35','#374527','#8BB400','#797868','#C6005A','#3B000A','#C86240','#29607C','#402334','#7D5A44','#CCB87C','#B88183','#AA5199','#B5D6C3','#A38469','#9F94F0','#A74571','#B894A6','#71BB8C','#00B433','#789EC9','#6D80BA','#953F00','#5EFF03','#E4FFFC','#1BE177','#BCB1E5','#76912F','#003109','#0060CD','#D20096','#895563','#29201D','#5B3213','#A76F42','#89412E','#1A3A2A','#494B5A','#A88C85','#F4ABAA','#A3F3AB','#00C6C8','#EA8B66','#958A9F','#BDC9D2','#9FA064','#BE4700','#658188','#83A485','#453C23','#47675D','#3A3F00','#061203','#DFFB71','#868E7E','#98D058','#6C8F7D','#D7BFC2','#3C3E6E','#D83D66','#2F5D9B','#6C5E46','#D25B88','#5B656C','#00B57F','#545C46','#866097','#365D25','#252F99','#00CCFF','#674E60','#FC009C','#92896B'}); 
% lab = rgb2lab(col);
% col = col(lab(:,1)>=5 & lab(:,1)<=95,:); % gets rid of very light or dark;  
col = mat2rgb(lon,cmocean('phase'),[-1 1]*180); 

col(183,:) = [0 0 0];

lb = [183 121 61 39 33 ];% include 60 for George Vi
lw1 = 0.2; 
lw2 = 0.4; 

vp = 0.035; 
hp = 0.02; 
xcol = .15*[1 1 1]; 
ycol = 0.6*[1 1 1]; 
figure('pos',[200 100 680 290])

%ax=subplot(2,4,1);
ax=subsubplot(1,2,1,2,2,1,'vpad',vp,'hpad',hp);
hold on
for k=1:183
   if ismember(k,lb)
      pl(k) = plot(D.year,glf_tp(k,:)-glf_tp(k,1),'color',col(k,:),'linewidth',lw2); 
      if k==183         
         txt(k)=text(D.year(13),glf_tp(k,13)-glf_tp(k,1),Dn.name{k},'fontsize',6,'fontweight','bold','horiz','right','vert','bot','color',col(k,:));
      else
         txt(k)=text(D.year(end),glf_tp(k,end)-glf_tp(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:));
      end
   else         
      pl(k) = plot(D.year,glf_tp(k,:)-glf_tp(k,1),'color',col(k,:),'linewidth',lw1); 
   end
end
pl(183).LineWidth = 1; 
%txt(183).FontWeight='bold'; 
%txt(48).VerticalAlignment='top'; 
txt(61).VerticalAlignment='bot'; 
box off
axis tight
title('Response to Thinning','fontsize',8,'fontweight','normal')
ylim([-2.5 19])
set(gca,'fontsize',7,'xtick',1995:5:2020,'xcolor',xcol,'color','none')
ntitle(' a','location','nw','fontweight','bold','fontsize',8,'horiz','right')
xtickangle(0)
set(gca,'ytick',-5:5:25,'yticklabel','')

tl = 5:5:10; 
text(D.year(1)+zeros(size(tl)),tl,num2str(tl'),'fontsize',7,'color',ycol,'vert','bot')
text(D.year(1),15,'15 Gt yr^{-1}','fontsize',7,'color',ycol,'vert','bot')
hl=hline(0:5:15,'color',ycol,'linewidth',.2); 
uistack(hl(:),'bottom');
set(gca,'ycolor','none')
ylabel('Change in GL flux','color',ycol,'visible','on')


lb=[183 121 153];
%ax(2) = subplot(2,4,2);
ax(2)=subsubplot(1,2,1,2,2,2,'vpad',vp,'hpad',hp);
hold on
for k=1:183
   if ismember(k,lb)
      pl(k) = plot(I.year,glf_cp(k,:)-glf_cp(k,1),'color',col(k,:),'linewidth',lw2); 

      if k==183         
         txt(k)=text(I.year(23),glf_cp(k,23)-glf_cp(k,1),Dn.name{k},'fontsize',6,'fontweight','bold','horiz','right','vert','bot','color',col(k,:));
      else
         txt(k)=text(I.year(end),glf_cp(k,end)-glf_cp(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:));
      end
   else
      pl(k) = plot(I.year,glf_cp(k,:)-glf_cp(k,1),'color',col(k,:),'linewidth',lw1); 
   end
end
pl(183).LineWidth = 1; 
txt(183).FontWeight='bold'; 
box off

title('Response to Calving','fontsize',8,'fontweight','normal')
ylim([-2.5 19])
set(gca,'fontsize',7,'xtick',1995:5:2020,'ycolor','none','xcolor',xcol,'color','none')
ntitle(' b','location','nw','fontweight','bold','fontsize',8,'horiz','right')
xtickangle(0)
hl=hline(0:5:15,'color',ycol,'linewidth',.2); 
uistack(hl(:),'bottom');

lb = [183 132 48 10];

tmp = glf_tf2; 
tmp(:,end) = glf_cf(:,end); 
%ax(3) = subplot(2,4,5);

ax(3)=subsubplot(1,2,1,2,2,3,'vpad',vp,'hpad',hp);

haha=1:1:101; 

hold on
for k=1:183
   %sc = (glf_cf(k,end)-glf_cf(k,1))/(tmp(k,end)-tmp(k,1)); 
   sc=1; 
   if ismember(k,lb)
      pl(k) = plot(haha-1,(tmp(k,haha)-tmp(k,1))*sc,'-','color',col(k,:),'linewidth',lw2); 
      if k==183
         txt(k)=text(80,(tmp(k,80)-tmp(k,1))*sc,Dn.name{k},'fontsize',6,'fontweight','bold','horiz','right','vert','bot','color',col(k,:));
      else
         txt(k)=text(100,(tmp(k,end)-tmp(k,1))*sc,Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:));
      end
   else 
      pl(k) = plot(haha-1,(tmp(k,haha)-tmp(k,1))*sc,'-','color',col(k,:),'linewidth',lw1); 
   end
end
pl(183).LineWidth = 1; 
txt(183).FontWeight = 'bold'; 
txt(48).VerticalAlignment='bottom'; 
box off
axis([0 100 0 2500])
set(gca,'fontsize',7,'xtick',0:20:100,'color','none')
xtickangle(0)
xlabel 'Percent mass lost to thinning'
ntitle(' c','location','nw','fontweight','bold','fontsize',8,'horiz','right')
set(gca,'ytick',0:500:3000,'yticklabel','')
tl = 500:500:1500; 
text(zeros(size(tl)),tl,num2str(tl'),'fontsize',7,'color',ycol,'vert','bot')
text(0,2000,'2000 Gt yr^{-1}','fontsize',7,'color',ycol,'vert','bot')
hl=hline(0:500:2000,'color',ycol,'linewidth',.2); 
uistack(hl(:),'bottom');
set(gca,'ycolor','none')
ylabel('Change in GL flux','color',ycol,'visible','on')

lb=[183.00        132.00         48.00         10.00 ];
%ax(4) = subplot(2,4,6);
ax(4)=subsubplot(1,2,1,2,2,4,'vpad',vp,'hpad',hp);
hold on
for k=1:183
   pct = 100*(HM.M_calving_future(1,k)-HM.M_calving_future(:,k))./HM.M_calving_future(1,k);
   if ismember(k,lb)
      pl(k) = plot(pct(haha),glf_cf(k,haha)-glf_cf(k,1),'-','color',col(k,:),'linewidth',lw2); 
      if k==183
         txt(k)=text(pct(99),glf_cf(k,99)-glf_cf(k,1),Dn.name{k},'fontsize',6,'fontweight','bold','horiz','right','vert','bot','color',col(k,:));
      else
         txt(k)=text(100,glf_cf(k,end)-glf_cf(k,1),Dn.name{k},'fontsize',6,'horiz','left','vert','middle','color',col(k,:));
      end
   else
      pl(k) = plot(pct(haha),glf_cf(k,haha)-glf_cf(k,1),'-','color',col(k,:),'linewidth',lw1); 
   end
end
pl(183).LineWidth = 1; 
%txt(183).FontWeight = 'bold'; 
txt(48).VerticalAlignment='bottom'; 
%txt(155).VerticalAlignment='top'; 
box off
axis([0 100 0 2500])
set(gca,'fontsize',7,'xtick',0:20:100,'ycolor','none','color','none')
xlabel 'Percent mass lost to calving'
xtickangle(0)
ntitle(' d','location','nw','fontweight','bold','fontsize',7,'horiz','right')
hl=hline(0:500:2000,'color',ycol,'linewidth',.2); 
uistack(hl(:),'bottom');

ax(5)=subplot(1,2,2); 
scatter(glf_cf(1:181,1),dg(1:181),Asc,lon,'filled')
axis tight
box off
caxis([-180 180]) 
cmocean phase 
ind = glf_cf(:,end)>25; 
ind(dg>200) = true; 
ind([19 33 49 74 81 87 113 125 146 147 182:183]) = false; % remove labels that clutter 
text(glf_cf(ind,1),dg(ind),Dn.name(ind),'horiz','center','vert','bot','fontsize',6,...
   'fontangle','italic','color',.35*[1 1 1])
axis tight
title('Response to Total Ice Shelf Collapse','fontsize',8,'fontweight','normal')
xlabel('Modeled present-day GL flux (Gt yr^{-1})','fontsize',7)
ylabel('Acceleration (%)','fontsize',7)
set(gca,'fontsize',7)

[X,Y] = meshgrid(0:10:260,0:10:400);

%G = (1+Y/100).*X; 
G = (Y/100).*X; 
hold on

[C,hC] = contour(X,Y,G,0:50:1500,'k'); 
hC.Color = ycol;
hC.LineWidth = 0.2; 

uistack(hC,'bottom')
axis([0 260 0 380])

n = 50:50:950; 
ny = 100*n./(260*ones(size(n))); 
str = num2str(n','%-4.f'); 
txt = text(260*ones(size(n)),ny,str,'vert','middle','fontsize',6,'color',ycol); 
%txt2 = text(250,20,'50 Gt yr^{-1}','vert','middle','fontsize',6,'color',.8*[1 1 1]); 
txt2 = text(270,mean(ylim),'Instantaneous increase in GL flux Gt yr^{-1}',...
   'rotation',90,'fontsize',6,'color',ycol,'vert','top','horiz','center');
col = mat2rgb(lon,cmocean('phase'),[-1 1]*180); 
pos = get(gca,'outerposition');
set(gca,'outerposition',pos+[-.04 0 0 0]) % snuggle up
%
ntitle(' e ','location','nw','fontweight','bold','fontsize',8)

gp = plotboxpos(gca);
%axes('position',[gp(1)+(.75)*gp(3) gp(2) .25*gp(3) .25*gp(4)])
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

% export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_results_compiled.jpg','-pdf','-r600','-painters')
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/fig_4.eps','contenttype','vector')
%% Old stuff 

%% Generate pdfs of response to hypothetical thinning
% IGNORE THIS SECTION
% We gave up on the "years of continued thinning" approach. 
% 
% for kk=0:6
%    %figure('pos',[40 40 560 760])
%    clf
%    for k = 1:27
%       if (kk*27+k)<=183
%          subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 
% 
%          set(gca,'fontsize',5) 
%             plot(future_years,glf_tf(kk*27+k,:),'color',colt,'linewidth',lw2,'markersize',ms)
%             hold on
%             %plot(0:100,glf_cf(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)
% 
%             box off
%             axis tight
%             %xlim([1997 2022])
%             ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')
% 
%             ax = gca; 
%             %ax.YAxis.Exponent = 0;
%             %ytickformat('%.0f');
%             set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])
% 
%       end
%       if kk==6
%          k=21; % just to set the legend on the last plot properly
%       end
%    end
% 
%    % Make a legend for the last axes: 
%    lg = legend('response to thinning','location','southwest');
% 
%    % Create a new axis just to get its position: 
%    tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);
% 
%    % Move the legend to the empty axis position and delete the empty axes: 
%    lg.Position=tmpax.Position;
%    delete(tmpax)
% 
%    sgtitle('Modeled grounding line flux (Gt yr^{-1}) response to years of continued trends','fontsize',8) 
% 
%    export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_hypothetical_thinning_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')
% 
% end
% 
% %% Generate pdfs of response to hypothetical calving
% 
% for kk=0:6
%    figure('pos',[40 40 560 760])
%    for k = 1:27
%       if (kk*27+k)<=183
%          subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 
% 
%          set(gca,'fontsize',5) 
%             plot(0:100,glf_cf(kk*27+k,:),'color',colc,'linewidth',lw2,'markersize',ms)
%             hold on
%             %
%             box off
%             axis tight
%             %xlim([1997 2022])
%             ntitle(Dn.name{kk*27+k},'fontsize',5,'color','k')
% 
%             ax = gca; 
%             %ax.YAxis.Exponent = 0;
%             %ytickformat('%.0f');
%             set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])
% 
%       end
%       if kk==6
%          k=21; % just to set the legend on the last plot properly
%       end
%    end
% 
%    % Make a legend for the last axes: 
%    lg = legend('response to calving','location','southwest');
% 
%    % Create a new axis just to get its position: 
%    tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);
% 
%    % Move the legend to the empty axis position and delete the empty axes: 
%    lg.Position=tmpax.Position;
%    delete(tmpax)
% 
%    sgtitle('Modeled grounding line flux (Gt yr^{-1}) response to ice shelf loss, by percent area','fontsize',8) 
% 
%    export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_hypothetical_calving_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')
% 
% end
% 
% %% Total ice shelf collapse scenario 
% 
% % Change in GL flux (percent): 
% dg = 100*(glf_cf(:,end)-glf_cf(:,1))./glf_cf(:,1);
% 
% Asc = sqrt(Dn.area_km2); 
% Asc = 90*Asc/max(Asc); 
% 
% [~,lon] = ps2ll(Dn.x_center,Dn.y_center);
% 
% figure('pos',[200 200 308 231])
% scatter(glf_cf(1:181,1),dg(1:181),Asc,lon,'filled')
% axis tight
% box off
% caxis([-180 180]) 
% cmocean phase 
% ind = glf_cf(:,end)>25; 
% ind(dg>200) = true; 
% ind([19 33 49 74 81 87 113 146 147 182:183]) = false; % remove labels that clutter 
% text(glf_cf(ind,1),dg(ind),Dn.name(ind),'horiz','center','vert','bot','fontsize',6,...
%    'fontangle','italic','color',.4*[1 1 1])
% axis tight
% title('Total ice shelf collapse scenario','fontsize',8)
% xlabel('Modeled present-day GL flux (Gt/yr)','fontsize',7)
% ylabel('Percent change in GL flux','fontsize',7)
% set(gca,'fontsize',7)
% 
% [X,Y] = meshgrid(0:10:250,0:10:450);
% 
% %G = (1+Y/100).*X; 
% G = (Y/100).*X; 
% hold on
% 
% [C,hC] = contour(X,Y,G,0:50:1500,'k'); 
% hC.Color = .8*[1 1 1];
% hC.LineWidth = 0.3; 
% 
% uistack(hC,'bottom')
% axis([0 250 0 450])
% 
% n = 50:50:1100; 
% ny = 100*n./(250*ones(size(n))); 
% str = num2str(n','%-4.f'); 
% txt = text(250*ones(size(n)),ny,str,'vert','middle','fontsize',6,'color',.8*[1 1 1]); 
% col = mat2rgb(lon,cmocean('phase'),[-1 1]*180); 
% 
% gp = plotboxpos(gca);
% axes('position',[gp(1)+(.75)*gp(3) gp(2) .25*gp(3) .25*gp(4)])
%  
% hold on
% for k = 1:181
%    plot(P(k),'facecolor',col(k,:),'facealpha',1,'edgecolor','none')
% end
% 
% hant = antbounds('coast','polyshape','facecolor','w','facealpha',1,'edgecolor','none');
% uistack(hant,'bottom'); 
% axis tight off
% bedmachine('coast','color',0.3*[1 1 1],'linewidth',0.1)
% bedmachine('gl','color',0.3*[1 1 1],'linewidth',0.1)
% gp2 = plotboxpos(gca); 
% set(gca,'xcolor','none','ycolor','none','pos',[gp(1)+(.75)*gp(3) gp(2) gp2(3) gp2(4)])
% set(gca,'fontsize',7)
% 
% % export_fig('/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/catastrophic_calving_response.jpg','-pdf','-r600','-painters')


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
