
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
[xc,yc] = polycenter(D.x,D.y);
   
warning('off')
for k = 1:181
   P(k) = polyshape(D.x{k},D.y{k}); 
end
A = area(P); 

xc(A<229092155) = nan; 

xc([28 88 91 165 166 65 69 78 59 116 98 25 68 23 24 15 46 51 117 154 148 179]) = nan; % names we don't care about

D.name{84} = 'Larsen A'; 
D.name{85} = 'Larsen B'; 
D.name{86} = 'Larsen C'; 
D.name{87} = 'Larsen D'; 

xc(48) = -790240; % manual Filchner 
yc(48) = 945473; 

%% Make a map: 

% Colors for each year's coastlie
col = mat2rgb(year,turbo(276)); 

shelf = {'larsens','fris','amery',...
   'wilkins','center','shackleton',...
   'ase','ross','oates'}; 

figure('pos',[10 200 670 516])

pad = 0.003;

for kk = [1 2 3 4 6 7 8 9]
  
   subsubplot(3,3,kk,'vpad',pad,'hpad',pad); 
   
   hold on
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
         axis([1924094.26    2289209.67     535365.14     831726.35])
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
         axis([-507558      507551    -1796747     -972790])
         %axis([-501528.81     277001.60   -1590125.01    -958197.22]) 
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
   
   q = itslive_quiver; 
   q.Color = hex2rgb('d1b057'); 
   q.LineWidth = 0.2; 
   uistack(q,'bottom')
   
   m(kk) = modismoaps('contrast','w','year',2004); 
   lt=text(xc,yc,D.name,'horiz','center','vert','middle',...
      'fontsize',5,'clipping','on','color',0.3*[1 1 1],'fontangle','italic'); 
   lt(121).VerticalAlignment='bottom'; %pine island
   lt(153).VerticalAlignment='bottom'; % thwaites
   lt(113).VerticalAlignment='bottom'; % ninnis
   lt(105).VerticalAlignment='bottom'; % mertz
   
   drawnow
   letters = {' a ',' b ',' c ',' d ',' e ',' f ',' g ',' h ',' i '}; 
   nt(kk) = ntitle(letters{kk},'location','nw','fontsize',8,'fontweight','bold',...
      'color',1*[1 1 1],'backgroundcolor','k','margin',1); 
   
end

ax(5) = subsubplot(3,3,5,'vpad',pad,'hpad',pad);
hold on
for k = 1:24
   h(k) = plot(cx{k},cy{k},'color',col(k,:),'linewidth',0.3); 
end
axis tight
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
uistack(hb,'bottom')
axis off 

for kk = [1 2 3 4 6 7 8 9]
   
   plot([xl(kk,1) xl(kk,2) xl(kk,2) xl(kk,1) xl(kk,1)],...
      [yl(kk,1) yl(kk,1) yl(kk,2) yl(kk,2) yl(kk,1)],...
      'color',.5*[1 1 1],'linewidth',0.2)
   
      switch letters{kk}
         case ' a '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',.5*[1 1 1],...
               'horiz','right','vert','top')
         case {' b ',' i ',' f '}
            text(xl(kk,1),yl(kk,2),letters{kk},'fontsize',5,'color',.5*[1 1 1],...
               'horiz','left','vert','top')
         case {' g ',' h '}
            text(xl(kk,1),yl(kk,1),letters{kk},'fontsize',5,'color',.5*[1 1 1],...
               'horiz','left','vert','bot')
         case ' d '
            text(xl(kk,1),yl(kk,1),letters{kk},'fontsize',5,'color',.5*[1 1 1],...
               'horiz','right','vert','bot')
         case ' c '
            text(xl(kk,2),yl(kk,2),letters{kk},'fontsize',5,'color',.5*[1 1 1],...
               'horiz','left','vert','top')
            
      end
end
axis tight


tcb = textcolorbar(unique(floor(year)),'colormap',turbo(256),'location','center','fontsize',5); 

tcb2 = text(825169+240e3,1394090,tcb.String(1:11,:),'fontsize',5,'vert','top','horiz','right');
tcb3 = text(825169+240e3,1394090,tcb.String(12:end,:),'fontsize',5,'vert','top','horiz','left');
delete(tcb)

m(5) = modismoaps('contrast','w','year',2004); 
set(gcf,'color','k') 

% export_fig test2.png -r600 -p0.005
%% * * * * * * SUBFUNCTIONS * * * * * * 

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

