
% This script creates a single-figure compilation of maps of Antarctic coastal 
% change 1997-2021. 
% 
% Chad Greene, October 2021. 

%% Load data

load('icemask_composite.mat','cx','cy','year')

% fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
% vx = permute(h5read(fn,'/vx'),[2 1]); 
% vy = permute(h5read(fn,'/vy'),[2 1]); 

%% Calculate center locations of ice shelf polygons. 
% Because centroids are outside the actual ice shelf outline some places,
% like crescent-shaped George VI. 

D = load('iceshelves_2008_v2.mat');
   
xc = D.x_center; 
yc = D.y_center; 

xc(D.area_km2<229.092155) = nan; 

xc([28 88 91 165 166 65 69 78 59 116 98 25 68 23 24 15 46 51 117 154 148 179 111 99 41 20 104]) = nan; % names we don't care about

D.name{84} = 'Larsen A'; 
D.name{85} = 'Larsen B'; 
D.name{86} = 'Larsen C'; 
D.name{87} = 'Larsen D'; 

xc(48) = -790240; % manual Filchner 
yc(48) = 945473; 

%% Make a map: 

cmp = 'tropical'; 
switch cmp 
   case 'pride'
      load('pride.mat')

      lch_pride = colorspace('rgb->lch',pride);

      tmp = linspace(18,95,256); 
      lch_pride(1:256,1) = tmp; 
      lch_pride(257:end,1) = tmp(255:-1:1); 

      mm = movmean(lch_pride(:,1),75); 
      lch_pride(256-100:256+100,1) = mm(256-100:256+100); 

      lch_pride(:,2) = movmean(lch_pride(:,2),75)+5;

      cm = colorspace('lch->rgb',lch_pride); 
      
   case 'neon'
      load neon
      cm = flipud(neon); 
      
   case 'tropical'
      load tropical
      cm = flipud(tropical);
      
   case 'pepper'
      load pepper
      cm = flipud(pepper); 
      
   case 'parula'
      lch_parula = colorspace('rgb->lch',parula(256));
      lch_parula(:,1) = linspace(20,90,256); 
      cm=flipud(colorspace('lch->rgb',lch_parula)); 
      
   case 'chroma'
      load chroma
      cm = flipud(chroma(50:230,:)); 
end

% Colors for each year's coastlie
col = mat2rgb(year,cm); 

shelf = {'larsens','fris','amery',...
   'wilkins','center','shackleton',...
   'ase','ross','oates'}; 

figure('pos',[10 200 670 516])

pad = 0.003;

for kk = [1 2 3 4 6 7 8 9]
  
   subsubplot(3,3,kk,'vpad',pad,'hpad',pad); 
   
   hold on
   % Plot an undercarraige: 
   if ismember(kk,[2 3 8 9])
   for k = 24
      h(k) = plot(cx{k},cy{k},'w','linewidth',0.6); 
   end
   end

   for k = 1:24
      h(k) = plot(cx{k},cy{k},'color',col(k,:),'linewidth',0.3); 
   end
   axis tight
   hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
   uistack(hb,'bottom')
   axis off 

   switch shelf{kk}
      
      case 'amery'
         %axis([1907706.00    2299893.00     543775.00     862110.00])
         %axis([1924094.26    2289209.67     535365.14     831726.35])
         axis([1946677.39    2285375.92     547598.37     822517.21])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'harald'
         axis([1196007.00    1552509.00    1678261.00    1967631.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'fimbul'
         axis([-371588.00     365766.00    1832896.00    2431401.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'brunt'
         axis([ -992943.00    -188319.00    1330067.00    1983176.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','length',100); 

      case 'fris'
         axis([ -1491530.00    -729114.00     580480.00    1199328.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','length',100); 

      case 'larsens'
         % axis([-2434759.00   -1808135.00     931420.00    1440047.00])
         axis([  -2444448.42   -1969057.53    1007510.13    1393382.11])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se'); 

      case 'wilkins'
         %axis([ -2228118    -1698659      389539      819297])
         axis([ -2211220.75   -1717553.16     441721.01     842427.38])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'abbot'
         axis([-2193626    -1661698     -464662      -32900])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'ase'
         axis([-1931795    -1420951     -699668     -285019])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'getz'
         axis([-1684886    -1099570    -1186787     -711690])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'sulzberger'
         axis([ -1059894     -513929    -1469287    -1026130])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'ross'
         %axis([-507558      507551    -1796747     -972790])
         axis([-506341.08     332780.77   -1548399.94    -867290.50])
         
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'oates'
         axis([1020867     1476174    -2251537    -1881967])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'porpoise'
         axis([1835885     2205088    -1733106    -1433426])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','color','w');

      case 'sabrina'
         axis([2034711     2451438    -1401418    -1063164])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','color','w'); 

      case 'vincennes'
         axis([2315562     2507597     -927267     -771394])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','length',25); 

      case 'shackleton'
         %axis([2299078     2867870     -672477     -210792])
         axis([2414909.15    2790235.81    -548289.40    -243639.04])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','k'); 

      case 'west'
         axis([2309813     2726932       44406      382979])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 
   end
   
   hsb(1).LineWidth = 1; 
   
   xl(kk,:) = xlim; 
   yl(kk,:) = ylim; 
   
   q(kk) = itslive_quiver; 
   q(kk).LineWidth = 0.2; 
   uistack(q(kk),'bottom')
   
   m(kk) = modismoaps('contrast','w','year',2004); 
%    lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
%       'fontsize',5,'clipping','on','color',hex2rgb('af8b69'),'fontangle','italic'); 
%    lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
%       'fontsize',5,'clipping','on','color',hex2rgb('4d7798'),'fontangle','italic'); 
   lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
      'fontsize',5,'clipping','on','color',hex2rgb('c26a4b'),'fontangle','italic');    
   
   lt(121).VerticalAlignment='bottom'; %pine island
   lt(153).VerticalAlignment='bottom'; % thwaites
   lt(153).HorizontalAlignment='left'; % thwaites
   lt(113).VerticalAlignment='bottom'; % ninnis
   lt(105).VerticalAlignment='bottom'; % mertz
   
   drawnow
   letters = {' a ',' b ',' c ',' d ',' e ',' f ',' g ',' h ',' i '}; 
   nt(kk) = ntitle(letters{kk},'location','nw','fontsize',8,'fontweight','bold',...
      'color',1*[1 1 1],'backgroundcolor','k','margin',1); 
   
end

ax(5) = subsubplot(3,3,5,'vpad',pad,'hpad',pad);
hold on

% % Plot an undercarraige: 
% for k = 1:24
%    h(k) = plot(cx{k},cy{k},'w','linewidth',0.8); 
% end

for k = 1:24
   h(k) = plot(cx{k},cy{k},'color',col(k,:),'linewidth',0.2); 
end
axis tight
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
uistack(hb,'bottom')
axis off 

%insetcol = hex2rgb('#9e7da5'); 
insetcol = hex2rgb('dd99c1'); 
for kk = [1 2 3 4 6 7 8 9]
   
   plot([xl(kk,1) xl(kk,2) xl(kk,2) xl(kk,1) xl(kk,1)],...
      [yl(kk,1) yl(kk,1) yl(kk,2) yl(kk,2) yl(kk,1)],...
      'color',insetcol,'linewidth',0.2)
   
      switch letters{kk}
         case ' a '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','top')
         case {' b ',' i ',' f '}
            text(xl(kk,1),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','top')
         case {' g ',' h '}
            text(xl(kk,1),yl(kk,1),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','bot')
         case ' d '
            text(xl(kk,1),yl(kk,1),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','bot')
         case ' c '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','top')
            
      end
end
axis tight

tcb = textcolorbar(unique(floor(year)),'colormap',cm,'location','center','fontsize',5); 

tcb2 = text(825169+190e3,1394090,tcb.String(1:11,:),'fontsize',5,'vert','top','horiz','right');
tcb3 = text(825169+190e3,1394090,tcb.String(12:end,:),'fontsize',5,'vert','top','horiz','left');
delete(tcb)


S = shaperead('/Users/cgreene/Downloads/FlowLines_InteriorSeeds_NoPH_FILLED/FlowLines_InteriorSeeds_NoPH_FILLED.shp');

q(5) = plot([S.X],[S.Y],'color','k','linewidth',0.2); 
% q = itslive_quiver('density',150); 
% q.Color = hex2rgb('77c1be'); 
% q.LineWidth = 0.2; 
% q.AutoScaleFactor=4; 
uistack(q(5),'bottom')
   
m(5) = modismoaps('contrast','w','year',2004); 
set(gcf,'color','k') 

% Flowline and vector color: 
for kk=1:9;q(kk).Color(1:3)=hex2rgb('69b1e2'); end
for kk=1:9;q(kk).Color(4) = 0.5; end

%export_fig coastal_change_maps.jpg -r600 -painters -p0.005
%export_fig coastal_change_maps_1200dpi.jpg -r1200 -painters -p0.005
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/fig_1.eps','contenttype','vector')
%% B-list maps 
% These are the second-rate or otherwise uninteresting ice shelves that
% didn't make the cut for the primary map. 


% Reload ice shelves without deleted bits: 
D = load('iceshelves_2008_v2.mat');
xc = D.x_center; 
yc = D.y_center; 
xc(D.area_km2<100) = nan; 

shelf = {'brunt','fimbul','harald',...
   'abbot','center','west',...
   'getz','drygalski','sabrina'}; 

figure('pos',[10 200 670 516])

pad = 0.003;

for kk = [1 2 3 4 6 7 8 9]
  
   subsubplot(3,3,kk,'vpad',pad,'hpad',pad); 
   
   hold on
   % Plot an undercarraige: 
   if ismember(kk,[1 6 7 8 9])
   for k = 24
      h(k) = plot(cx{k},cy{k},'w','linewidth',0.6); 
   end
   end

   for k = 1:24
      h(k) = plot(cx{k},cy{k},'color',col(k,:),'linewidth',0.3); 
   end
   axis tight
   hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
   uistack(hb,'bottom')
   axis off 

   switch shelf{kk}
      
      case 'amery'
         %axis([1907706.00    2299893.00     543775.00     862110.00])
         %axis([1924094.26    2289209.67     535365.14     831726.35])
         axis([1946677.39    2285375.92     547598.37     822517.21])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'harald'
         axis([1196007.00    1552509.00    1678261.00    1967631.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'drygalski'
         %axis([235733.15     536781.65   -1717608.99   -1473249.79])
         axis([361645.68     535685.69   -1631342.88   -1490075.69])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 
         
      case 'fimbul'
         axis([-371588.00     365766.00    1832896.00    2431401.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 

      case 'brunt'
         axis([ -992943.00    -188319.00    1330067.00    1983176.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','length',100); 

      case 'fris'
         axis([ -1491530.00    -729114.00     580480.00    1199328.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','length',100); 

      case 'larsens'
         % axis([-2434759.00   -1808135.00     931420.00    1440047.00])
         axis([  -2444448.42   -1969057.53    1007510.13    1393382.11])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se'); 

      case 'wilkins'
         %axis([ -2228118    -1698659      389539      819297])
         axis([ -2211220.75   -1717553.16     441721.01     842427.38])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'abbot'
         axis([-2193626    -1661698     -464662      -32900])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'ase'
         axis([-1931795    -1420951     -699668     -285019])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'getz'
         axis([-1684886    -1099570    -1186787     -711690])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'sulzberger'
         axis([ -1059894     -513929    -1469287    -1026130])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'ross'
         %axis([-507558      507551    -1796747     -972790])
         axis([-506341.08     332780.77   -1548399.94    -867290.50])
         
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'oates'
         axis([1020867     1476174    -2251537    -1881967])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','w'); 

      case 'porpoise'
         axis([1835885     2205088    -1733106    -1433426])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','color','w');

      case 'sabrina'
         axis([2034711     2451438    -1401418    -1063164])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','se','color','w'); 

      case 'vincennes'
         axis([2315562     2507597     -927267     -771394])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','length',25); 

      case 'shackleton'
         %axis([2299078     2867870     -672477     -210792])
         axis([2414909.15    2790235.81    -548289.40    -243639.04])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw','color','k'); 

      case 'west'
         axis([2309813     2726932       44406      382979])
         [hsb(1),hsb(2)] = scalebarps('fontsize',5,'location','sw'); 
   end
   
   hsb(1).LineWidth = 1; 
   
   xl(kk,:) = xlim; 
   yl(kk,:) = ylim; 
   
   q(kk) = itslive_quiver; 
   q(kk).LineWidth = 0.2; 
   uistack(q(kk),'bottom')
   
   m(kk) = modismoaps('contrast','w','year',2004); 
%    lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
%       'fontsize',5,'clipping','on','color',hex2rgb('af8b69'),'fontangle','italic'); 
%    lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
%       'fontsize',5,'clipping','on','color',hex2rgb('4d7798'),'fontangle','italic'); 
   lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
      'fontsize',5,'clipping','on','color',hex2rgb('c26a4b'),'fontangle','italic');    
   
%    lt(121).VerticalAlignment='bottom'; %pine island
%    lt(153).VerticalAlignment='bottom'; % thwaites
%    lt(153).HorizontalAlignment='left'; % thwaites
%    lt(113).VerticalAlignment='bottom'; % ninnis
%    lt(105).VerticalAlignment='bottom'; % mertz
   
   drawnow
   letters = {' a ',' b ',' c ',' d ',' e ',' f ',' g ',' h ',' i '}; 
   nt(kk) = ntitle(letters{kk},'location','nw','fontsize',8,'fontweight','bold',...
      'color',1*[1 1 1],'backgroundcolor','k','margin',1); 
   
end

ax(5) = subsubplot(3,3,5,'vpad',pad,'hpad',pad);
hold on

% % Plot an undercarraige: 
% for k = 1:24
%    h(k) = plot(cx{k},cy{k},'w','linewidth',0.8); 
% end

for k = 1:24
   h(k) = plot(cx{k},cy{k},'color',col(k,:),'linewidth',0.2); 
end
axis tight
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
uistack(hb,'bottom')
axis off 

%insetcol = hex2rgb('#9e7da5'); 
insetcol = hex2rgb('dd99c1'); 
for kk = [1 2 3 4 6 7 8 9]
   
   plot([xl(kk,1) xl(kk,2) xl(kk,2) xl(kk,1) xl(kk,1)],...
      [yl(kk,1) yl(kk,1) yl(kk,2) yl(kk,2) yl(kk,1)],...
      'color',insetcol,'linewidth',0.2)
   
      switch letters{kk}
         case {' a ',' d '}
            text(xl(kk,1),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','top')
         case ' c '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','top')
         case {' b '}
            text(xl(kk,1),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','top')
         case {' g '}
            text(xl(kk,1),yl(kk,1),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','left','vert','bot')
         case ' h '
            text(xl(kk,1),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','top')
         case ' i '
            text(xl(kk,2),yl(kk,1),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','bot')
         case ' f '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',insetcol,...
               'horiz','right','vert','top')
      end
end
axis tight

tcb = textcolorbar(unique(floor(year)),'colormap',cm,'location','center','fontsize',5); 

tcb2 = text(825169+190e3,1394090,tcb.String(1:11,:),'fontsize',5,'vert','top','horiz','right');
tcb3 = text(825169+190e3,1394090,tcb.String(12:end,:),'fontsize',5,'vert','top','horiz','left');
delete(tcb)


S = shaperead('/Users/cgreene/Downloads/FlowLines_InteriorSeeds_NoPH_FILLED/FlowLines_InteriorSeeds_NoPH_FILLED.shp');

q(5) = plot([S.X],[S.Y],'color','k','linewidth',0.2); 
% q = itslive_quiver('density',150); 
% q.Color = hex2rgb('77c1be'); 
% q.LineWidth = 0.2; 
% q.AutoScaleFactor=4; 
uistack(q(5),'bottom')
   
m(5) = modismoaps('contrast','w','year',2004); 
set(gcf,'color','k') 

% Flowline and vector color: 
for kk=1:9;q(kk).Color(1:3)=hex2rgb('69b1e2'); end
for kk=1:9;q(kk).Color(4) = 0.5; end

%export_fig coastal_change_maps_part2.jpg -r600 -painters -p0.005
%export_fig coastal_change_maps_part2_1200dpi.jpg -r1200 -painters -p0.005
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/fig_ED1.eps','contenttype','vector')

%% * * * * * * SUBFUNCTIONS * * * * * * 

function RGB = mat2rgb(val,cmap,limits)
% 
% 
% Chad Greene wrote this, July ç2016. 

if nargin==2
   limits = [min(val) max(val)]; 
end

gray = mat2gray(val,limits); 
ind = gray2ind(gray,size(cmap,1));
RGB = cmap(ind+1,:); 

end

