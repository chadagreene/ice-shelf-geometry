
load icemask_composite

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-09-02.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
H_source = permute(h5read(fn,'/thickness_source'),[2 1]); 

x = h5read(fn,'/x');
y = h5read(fn,'/y'); 

D = load('iceshelves_2008_v2.mat');
%%

always_ice = all(ice,3); 
never_ice = ~any(ice,3); 

%%

mask2 = cube2rect(mask,~always_ice & ~never_ice); 
ice2 = cube2rect(ice,~always_ice & ~never_ice); 
H2 = cube2rect(H,~always_ice & ~never_ice); 

A = NaN(24,181); % Not total area! Just the active  area.
M = A; 
for k = 1:181
   
   A(:,k) = sum(ice2(:,mask2==k),2)*240^2; 
   M(:,k) = sum(H2(:,mask2==k).*ice2(:,mask2==k),2)*917*1e-12*240^2;
   
end

Atr = trend(A,year)/1000^2;
Mtr = trend(M,year);

%%

figure
bar(A(end,:)-A(1,:))
axis tight
box off
set(gca,'xtick',1:181,'xticklabel',D.name)
ylabel 'trend km^2/yr'

figure
plot(A(end,:)-A(1,:),M(end,:)-M(1,:),'o')
hold on
text(A(end,:)-A(1,:),M(end,:)-M(1,:),D.name,'horiz','center')
xlabel 'area change (km^2)'
ylabel 'mass chage (Gt)'

%%

[X,Y] = meshgrid(x,y); 
tidal = istidal(X,Y); 

Atr_map = nan(size(mask)); 
Adiff_map = nan(size(mask)); 
Amax = nan(181,1); 
for k = 1:181
   Atr_map(tidal & mask==k & ~never_ice) = Atr(k); 
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