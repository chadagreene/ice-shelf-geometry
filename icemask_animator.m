
% This script creates maps and animations of Antarctic coastal 
% change 1997-2021. 
% 
% Chad Greene, October 2021. 

%% Load data

load icemask_composite

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-15.h5'; 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

%% Calculate center locations of ice shelf polygons. 
% Because centroids are outside the actual ice shelf outline some places,
% like crescent-shaped George VI. 

D = load('iceshelves_2008_v2.mat');
xc = nan(181,1); 
yc = xc; 
for k = 1:181
   P = polyshape(D.x{k},D.y{k},'Simplify',true); 
   
   % Get Delaunay: 
   T = triangulation(P);
   
   % Circumcenters are the dots at the center of the triangulations: 
   [C,r] = circumcenter(T);
   
   % Index of centerpoint with largest radius from boundary:  
   [~,ind] = max(r); 
   xc(k) = C(ind,1); 
   yc(k) = C(ind,2); 
end

%% Make a map: 

% Colors for each year's coastlie
col = mat2rgb(year,parula); 

figure
hold on
for k = 1:24
   h(k) = plot(cx{k},cy{k},'color',col(k,:)); 
end
axis tight
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
uistack(hb,'bottom')

axis off 
set(gca,'pos',[0 0 1 1])
cb = colorbar('location','west');
tmp = cb.Position; 
cb.Color = 0.5*[1 1 1];
tmp(2) = 0.03; 
tmp(3) = tmp(3)/2; 
tmp(4) = tmp(4)/2;

posll = tmp; 
posul = tmp; 
posul(2) = .5; 
posur = posul; 
posur(1) = 1-posul(1); 
poslr = posur; 
poslr(2) = posll(2); 

cb.Position = tmp;
cb.FontSize = 8;
set(cb,'ytick',[1997:4:2021])
caxis([1997 2021])

% export_fig /Users/cgreene/Documents/MATLAB/DEM_generation/coastal_change_maps/coastalchange_map_antarctica.png -r600 -p0.01 -painters

txt = text(xc,yc,D.name,'horiz','center','vert','middle','fontsize',7,'color',.2*[1 1 1],'fontangle','italic','clipping','on');

%% Save each regional map: 

regions = {'amery','harald','fimbul','brunt','fris','larsens','wilkins','abbot',...
   'ase','getz','sulzberger','ross','oates','porpoise','sabrina','vincennes','shackleton','west'}; 

for k = 1:length(regions)
   try 
      delete(hsb(:))
   end
   try
      delete(m)
   end

   switch regions{k}
      case 'amery'
         axis([1907706.00    2299893.00     543775.00     862110.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 
         cb.Position = posul; 

      case 'harald'
         axis([1196007.00    1552509.00    1678261.00    1967631.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 
         cb.Position = posur; 

      case 'fimbul'
         axis([-371588.00     365766.00    1832896.00    2431401.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 
         cb.Position = poslr; 

      case 'brunt'
         axis([ -992943.00    -188319.00    1330067.00    1983176.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','length',100); 
         cb.Position = posul; 

      case 'fris'
         axis([ -1491530.00    -729114.00     580480.00    1199328.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',100); 
         cb.Position = posul; 

      case 'larsens'
         axis([-2434759.00   -1808135.00     931420.00    1440047.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se'); 
         cb.Position = posur; 

      case 'wilkins'
         axis([ -2228118    -1698659      389539      819297])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posur; 

      case 'abbot'
         axis([-2193626    -1661698     -464662      -32900])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posul; 

      case 'ase'
         axis([-1931795    -1420951     -699668     -285019])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posur; 

      case 'getz'
         axis([-1684886    -1099570    -1186787     -711690])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posur; 

      case 'sulzberger'
         axis([ -1059894     -513929    -1469287    -1026130])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posur; 

      case 'ross'
         axis([-507558      507551    -1796747     -972790])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posur; 

      case 'oates'
         axis([1020867     1476174    -2251537    -1881967])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 
         cb.Position = posul; 

      case 'porpoise'
         axis([1835885     2205088    -1733106    -1433426])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w');
         cb.Position = posul;  

      case 'sabrina'
         axis([2034711     2451438    -1401418    -1063164])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w'); 
         cb.Position = posul; 

      case 'vincennes'
         axis([2315562     2507597     -927267     -771394])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',25); 
         cb.Position = poslr; 

      case 'shackleton'
         axis([2299078     2867870     -672477     -210792])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','k'); 
         cb.Position = poslr; 

      case 'west'
         axis([2309813     2726932       44406      382979])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 
         cb.Position = poslr; 

      otherwise
         axis tight
         error 'unrecognized ice shelf'

   end

   m = modismoaps('contrast','w','year',2004); 

   export_fig(['/Users/cgreene/Documents/MATLAB/DEM_generation/coastal_change_maps/coastalchange_map_',regions{k},'.png'],'-r600','-p0.01','-painters')

end

%% Calculate passive tracer positions

Ns = 2e6; % number of scattered points

Xs = NaN(Ns,length(year)); 
Ys = Xs; 

% These points are randomly spread across the whole map and more: 
Xs(:,1) = 2*2816512*(rand(Ns,1)-.5); 
Ys(:,1) = 2*2816512*(rand(Ns,1)-.5); 

for k = 2:length(year)
   dt = year(k)-year(k-1); 
   Xs(:,k) = Xs(:,k-1) + dt*interp2(x,y,vx,Xs(:,k-1),Ys(:,k-1));
   Ys(:,k) = Ys(:,k-1) + dt*interp2(x,y,vy,Xs(:,k-1),Ys(:,k-1));
end

for k = 1:length(year)
   ind = interp2(x,y,single(ice(:,:,k)),Xs(:,k),Ys(:,k),'nearest',0)==1; 
   Xs(~ind,k) = nan; 
   Ys(~ind,k) = nan; 
end

%% Animation 

figure('pos',[10 400 580 552])

hh = imagescn(x,y,ice(:,:,1));
hold on
axis tight off
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
ttl = ntitle(num2str(floor(year(1))),'fontweight','bold','color',.5*[1 1 1],'FontSize',12); 
col = cmocean('ice'); 
set(gcf,'color',col(1,:))
cmocean ice
txt = text(xc,yc,D.name,'horiz','center','vert','middle','fontsize',7,'color',.4*[1 1 1],'fontangle','italic','clipping','on');

hm = plot(Xs(:,1),Ys(:,1),'.','color',col(230,:),'markersize',1); 
set(gca,'pos',[0 0 1 1])
uistack(hh,'bottom')
uistack(hb,'top')
uistack(txt)
uistack(ttl,'top') 

col2 = mat2rgb(year,parula); 

for kk = 2:length(regions)
   try 
      delete(hsb(:))
   end
   try
      delete(h(:))
   end

   switch regions{kk}
      case 'amery'
         axis([1907706.00    2299893.00     543775.00     862110.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'harald'
         axis([1196007.00    1552509.00    1678261.00    1967631.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'fimbul'
         axis([-371588.00     365766.00    1832896.00    2431401.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'brunt'
         axis([ -992943.00    -188319.00    1330067.00    1983176.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','length',100); 

      case 'fris'
         axis([ -1491530.00    -729114.00     580480.00    1199328.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',100); 

      case 'larsens'
         axis([-2434759.00   -1808135.00     931420.00    1440047.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se'); 

      case 'wilkins'
         axis([ -2228118    -1698659      389539      819297])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'abbot'
         axis([-2193626    -1661698     -464662      -32900])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'ase'
         axis([-1931795    -1420951     -699668     -285019])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'getz'
         axis([-1684886    -1099570    -1186787     -711690])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'sulzberger'
         axis([ -1059894     -513929    -1469287    -1026130])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'ross'
         axis([-507558      507551    -1796747     -972790])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'oates'
         axis([1020867     1476174    -2251537    -1881967])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'porpoise'
         axis([1835885     2205088    -1733106    -1433426])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w');

      case 'sabrina'
         axis([2034711     2451438    -1401418    -1063164])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w'); 

      case 'vincennes'
         axis([2315562     2507597     -927267     -771394])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',25); 

      case 'shackleton'
         axis([2299078     2867870     -672477     -210792])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','k'); 

      case 'west'
         axis([2309813     2726932       44406      382979])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      otherwise
         axis tight
         error 'unrecognized ice shelf'
   end
   
   v = VideoWriter(['/Users/cgreene/Documents/MATLAB/DEM_generation/coastal_change_maps/coastalchange_map2_',regions{kk},'.mp4'],'MPEG-4');
   v.FrameRate = 6; 
   v.Quality = 100; 
   open(v)

   for k = 1:24
      hh.CData = ice(:,:,k); 
      ttl.String = num2str(floor(year(k))); 
      ttl.Color = col2(k,:); 
      hm.XData = Xs(:,k); 
      hm.YData = Ys(:,k); 
      
      % Make the last one semitransprent: 
      if k>1
         h(k-1).Color(4) = 0.5; 
      end
      % Plot the new one:
      h(k) = plot(cx{k},cy{k},'color',col2(k,:)); 
      
      drawnow
      frame = export_fig('-nocrop','-r500','-painters'); 
      writeVideo(v,frame)
      
   end

   close(v)
end

%% Animation without coastlines: 

figure('pos',[10 400 580 552])

hh = imagescn(x,y,ice(:,:,1));
hold on
axis tight off
hb = bedmachine('gl','color',.5*[1 1 1],'LineWidth',0.25); 
ttl = ntitle(num2str(floor(year(1))),'fontweight','bold','color',.5*[1 1 1],'FontSize',12); 
col = cmocean('ice'); 
set(gcf,'color',col(1,:))
cmocean ice
txt = text(xc,yc,D.name,'horiz','center','vert','middle','fontsize',7,'color',.4*[1 1 1],'fontangle','italic','clipping','on');

hm = plot(Xs(:,1),Ys(:,1),'.','color',col(230,:),'markersize',1); 
set(gca,'pos',[0 0 1 1])
uistack(hh,'bottom')
uistack(hb,'top')
uistack(txt)
uistack(ttl,'top') 

for kk = 1:length(regions)
   try 
      delete(hsb(:))
   end

   switch regions{kk}
      case 'amery'
         axis([1907706.00    2299893.00     543775.00     862110.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'harald'
         axis([1196007.00    1552509.00    1678261.00    1967631.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'fimbul'
         axis([-371588.00     365766.00    1832896.00    2431401.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      case 'brunt'
         axis([ -992943.00    -188319.00    1330067.00    1983176.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','length',100); 

      case 'fris'
         axis([ -1491530.00    -729114.00     580480.00    1199328.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',100); 

      case 'larsens'
         axis([-2434759.00   -1808135.00     931420.00    1440047.00])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se'); 

      case 'wilkins'
         axis([ -2228118    -1698659      389539      819297])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'abbot'
         axis([-2193626    -1661698     -464662      -32900])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'ase'
         axis([-1931795    -1420951     -699668     -285019])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'getz'
         axis([-1684886    -1099570    -1186787     -711690])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'sulzberger'
         axis([ -1059894     -513929    -1469287    -1026130])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'ross'
         axis([-507558      507551    -1796747     -972790])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'oates'
         axis([1020867     1476174    -2251537    -1881967])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','w'); 

      case 'porpoise'
         axis([1835885     2205088    -1733106    -1433426])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w');

      case 'sabrina'
         axis([2034711     2451438    -1401418    -1063164])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','se','color','w'); 

      case 'vincennes'
         axis([2315562     2507597     -927267     -771394])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','length',25); 

      case 'shackleton'
         axis([2299078     2867870     -672477     -210792])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw','color','k'); 

      case 'west'
         axis([2309813     2726932       44406      382979])
         [hsb(1),hsb(2)] = scalebarps('fontsize',8,'location','sw'); 

      otherwise
         axis tight
         error 'unrecognized ice shelf'

   end
   
   
   v = VideoWriter(['/Users/cgreene/Documents/MATLAB/DEM_generation/coastal_change_maps/coastalchange_map_',regions{kk},'.mp4'],'MPEG-4');
   v.FrameRate = 6; 
   v.Quality = 100; 
   open(v)


   for k = 1:24
      hh.CData = ice(:,:,k); 
      ttl.String = num2str(floor(year(k))); 
      hm.XData = Xs(:,k); 
      hm.YData = Ys(:,k); 
      drawnow
      frame = export_fig('-nocrop','-r500','-painters'); 
      writeVideo(v,frame)
      
   end

   close(v)
end

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

