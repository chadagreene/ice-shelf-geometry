


load icemask_composite

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-09-02.h5'; 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 

D = load('iceshelves_2008_v2.mat');

xc = nan(181,1); 
yc = xc; 
for k = 1:181
   P = polyshape(D.x{k},D.y{k}); 
   [xc(k),yc(k)] = centroid(P); 
end
%% 

col = mat2rgb(year,parula); 

figure
hold on
for k = 1:24
   h(k) = plot(cx{k},cy{k},'color',col(k,:)); 
end
axis tight
hb = bedmachine('gl','color',.5*[1 1 1]); 
uistack(hb,'bottom')

axis off 
cb = colorbar;
cb.FontSize = 8;
caxis([1997 2021])

txt = text(xc,yc,D.name,'horiz','center','vert','middle','fontsize',7,'color',.2*[1 1 1],'fontangle','italic','clipping','on');
uistack(txt,'bottom')
%%
try 
   delete(hsb(:))
end

switch region 
   case 'amery'
      axis([1907706.00    2299893.00     543775.00     862110.00])
      [hsb(1),hsb(2)] = scalebarps('location','sw'); 
      
   case 'harald'
      axis([1196007.00    1552509.00    1678261.00    1967631.00])
      [hsb(1),hsb(2)] = scalebarps('location','sw'); 
      
   case 'fimbul'
      axis([-371588.00     365766.00    1832896.00    2431401.00])
      [hsb(1),hsb(2)] = scalebarps('location','sw'); 
      
   case 'brunt'
      axis([ -992943.00    -188319.00    1330067.00    1983176.00])
      [hsb(1),hsb(2)] = scalebarps('location','se'); 
      
   case 'fris'
      axis([ -1491530.00    -729114.00     580480.00    1199328.00])
      [hsb(1),hsb(2)] = scalebarps('location','sw'); 
      
   case 'larsens'
      axis([-2434759.00   -1808135.00     931420.00    1440047.00])
      [hsb(1),hsb(2)] = scalebarps('location','se'); 
      
   otherwise
      axis tight
      
end
      
      

%%

Ns = 2e6; % number of scattered points

Xs = NaN(Ns,length(year)); 
Ys = Xs; 

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
   k
end

%%



%% 
k=1;
figure
hh = imagescn(x,y,ice(:,:,k)); 
hb=bedmachine('gl','color',[.5 .5 .5]);
cmocean ice 
hold on
col = cmocean('ice'); 
set(gcf,'color',col(1,:))
txt = ntitle(num2str(floor(year(k))),'fontweight','bold','color',.5*[1 1 1]); 
axis off 

hm = plot(Xs(:,1),Ys(:,1),'.','color',col(230,:),'markersize',1); 

uistack(hb,'top')
uistack(txt,'top') 
%%

gif('/Users/cgreene/Documents/MATLAB/DEM_generation/calving_amery.gif','DelayTime',1/6,'Resolution',300)
% 
% % Play the last frame 5 more times: 
for k = 1:5
   gif
end

for k = 2:24
   hh.CData = ice(:,:,k); 
   txt.String = num2str(floor(year(k))); 
   hm.XData = Xs(:,k); 
   hm.YData = Ys(:,k); 
   drawnow
   gif
end
% 
% % Play the last frame 5 more times: 
for k = 1:5
   gif
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

