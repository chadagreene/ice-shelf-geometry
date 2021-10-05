
load icemask_composite

% fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/itslive_coastmask_2021-10-01.h5'; 
% ice = permute(logical(h5read(fn,'/ice')),[2 1 3]); 
% year = h5read(fn,'/year'); 

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-09-02.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
H_source = permute(h5read(fn,'/thickness_source'),[2 1]); 


x = h5read(fn,'/x');
y = h5read(fn,'/y'); 
[X,Y] = meshgrid(x,y); 
tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]); 

D = load('iceshelves_2008_v2.mat');

%%
always_ice = all(ice,3); 
ever_ice = any(ice,3); 

%%

mask2 = cube2rect(mask,tidal & ever_ice); 
ice2 = cube2rect(ice,tidal & ever_ice); 
H2 = cube2rect(H,tidal & ever_ice); 

A = NaN(24,181); 
M = A; 
for k = 1:181
   
   A(:,k) = sum(ice2(:,mask2==k),2)*240^2; 
   M(:,k) = sum(H2(:,mask2==k).*ice2(:,mask2==k),2)*917*1e-12*240^2;
   
end

Adiff = A(end,:)-A(1,:);
CMdiff = M(end,:)-M(1,:);

%%

figure
bar(Adiff)
axis tight
box off
set(gca,'xtick',1:181,'xticklabel',D.name)
ylabel 'trend km^2/yr'

figure
plot(Adiff,Mdiff,'o')
hold on
text(Adiff,Mdiff,D.name,'horiz','center')
xlabel 'area change (km^2)'
ylabel 'mass chage (Gt)'

%%

yrf = h5read('FULL_CUBE_v4.h5','/t'); 
xf = h5read('FULL_CUBE_v4.h5','/x');
yf = h5read('FULL_CUBE_v4.h5','/y');
[Xf,Yf] = meshgrid(xf,yf); 

Hf = permute(h5read('FULL_CUBE_v4.h5','/H_filt10'),[3 2 1]); 

ind = yrf>1997; 
Hf = Hf(:,:,ind); 
yrf = yrf(ind); 
Htr = trend(Hf,yrf); 
Hdiff = Hf(:,:,end) - Hf(:,:,1); 

maskf = interp2(x,y,mask,Xf,Yf,'nearest'); 
always_ice_f = interp2(x,y,always_ice,Xf,Yf,'nearest',0); 
isf = isfinite(Htr);
HMdiff = nan(1,181); 
for k = 1:181
   HMdiff(k) = sum(Hdiff(isf & always_ice_f & maskf==k),'all')*917*1e-12*3000^2;
end

%%

figure
hp = plot(HMdiff,CMdiff,'.','color',rgb('black'));
txt(1) = text(-1600,0,{'T';'H';'I';'N';'N';'I';'N';'G'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(2) = text(1600,0,{'T';'H';'I';'C';'K';'E';'N';'I';'N';'G'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(3) = text(0,-1600,{'RETREAT'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(4) = text(0,1600,{'ADVANCE'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
uistack(hp,'top'); 

ind = abs(HMdiff)>50 | abs(CMdiff)>50;
text(HMdiff(ind),CMdiff(ind),D.name(ind),'horiz','center',...
   'vert','bottom','fontsize',7,'color',rgb('black'))
box off
axis tight
axis equal
axis([-1 1 -1 1]*1.01*max(abs(axis)))
xlabel('due to thinning or thickening (Gt)','fontsize',8)
ylabel('due to retreat or advance (Gt)','fontsize',8)
title 'total ice shelf mass change since 1997' 

hold on
p(1)=hline(0,'color',rgb('light gray'));
p(2)=vline(0,'color',rgb('light gray'));
%p(3) = plot(xlim,xlim,'color',rgb('light gray')); 
uistack(p,'bottom')
set(gca,'fontsize',8)

ax = axis; 
box on
set(gca,'color',rgb('pastel pink')); 

hp = patch([ax(1) ax(2) 0 ax(1)],[ax(3) ax(3) 0 ax(3)],'b'); 
hp.FaceColor = rgb('pastel blue'); 
hp.EdgeColor = 'none'; 
uistack(hp,'bottom') 

hp(2) = patch([ax(1) 0 ax(2) ax(1)],[ax(4) 0 ax(4) ax(4)],'b'); 
hp(2).FaceColor = rgb('pastel blue'); 
hp(2).EdgeColor = 'none'; 
uistack(hp(2),'bottom') 

% export_fig calving_vs_thinning_masschange.png -r600 -p0.01 

%% 

figure
plot(100*HMdiff./M(1,:),100*CMdiff./M(1,:),'.','color',rgb('dark red'))
text(100*HMdiff./M(1,:),100*CMdiff./M(1,:),D.name,'horiz','center',...
   'vert','bottom','fontsize',7,'color',rgb('dark red'),'clipping','on')
box off
axis tight
axis equal
axis([-1 1 -1 1]*100)
xlabel('due to thinning','fontsize',8)
ylabel('due to change in extents','fontsize',8)
title 'Percent change in mass since 1997'

hold on
p(1) = plot(xlim,xlim,'color',rgb('light gray')); 
p(2)=hline(0,'color',rgb('light gray'));
p(3)=vline(0,'color',rgb('light gray'));
uistack(p,'bottom')
set(gca,'fontsize',8)
%ylim([-80 50])


% export_fig calving_vs_thinning_masschange_percent.png -r600 -p0.01 


%%

CM_percent_map = nan(size(mask)); 

for k = 1:181
   CM_percent_map(mask==k & tidal & ever_ice) = 100*CMdiff(k)/(CMdiff(k) + HMdiff(k)); 
end

figure
imagescn(x,y,CM_percent_map)
caxis([0 100])
bedmachine('gl','color',[.5 .5 .5])
axis off
cmocean -bal

%%



Atr_map = nan(size(mask)); 
Adiff_map = nan(size(mask)); 
Amax = nan(181,1); 
for k = 1:181
   Atr_map(tidal & mask==k & ever_ice) = Atr(k); 
   Amax(k) = sum(tidal & mask==k & ~never_ice,'all')*240^2;
end
for k = 1:181
   
   Adiff_map(tidal & mask==k & ~never_ice) = A(end,k)-A(1,k); 
end

A_map = nan(size(mask)); 
for k = 1:181
   
   A_map(tidal & mask==k & ice(:,:,1)) = Amax(k); 
end

figure
h = imagescn(x,y,100*Adiff_map./A_map); 
colorbar
caxis([-1 1]*100)
cmocean -bal
bedmachine('gl','color',[.5 .5 .5])
axis off

%%
iceshelf = 1; 

msk = mask==iceshelf & ~always_ice & ~never_ice; 

Ai = local(ice,msk,@sum)*.24^2; 

figure
plot(R.yr+.75,AR - Ai(1),'o','linewidth',1)
hold on
plot(MOA.yr+.2,AMOA - Ai(1),'^','linewidth',1)
plot(S.yr+.2,AS - Ai(1),'p','linewidth',1)
plot(F.years+.2,AF - Ai(1),'s','linewidth',1)
plot(year,Ai - Ai(1),'k-','linewidth',2)

box off 
axis tight
set(gcf,'renderer','painters') 
ylabel 'total ice area change (km^2)' 
legend('Radarsat','MOA','Fraser Sentinel 1a','Fraser MODIS','composite','location','best')
legend boxoff 
title(D.name{iceshelf})

