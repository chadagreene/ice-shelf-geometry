

load('thickness_and_melt_rates.mat','MR','MR_err','SMB','SMB_err','dHdt','dHdt_err','readme')
load('calving_flux_timeseries.mat')
D = load('iceshelves_2008_v2.mat');

%% Add the 2021 Ronne calving to the time series: 

% Convert to thousands of square km: 
A_calving = A_calving*1e-9; 

ronne_mass_per_area = diff(A_calving(:,132))\diff(M_calving(:,132));

year = [year 2021.4]; 
A_calving = [A_calving;A_calving(end,:)];
A_calving(end,132) = A_calving(end-1,132)-4.320; % Iceberg A-76 https://www.esa.int/ESA_Multimedia/Images/2021/05/Meet_the_world_s_largest_iceberg#.YKTfep3wFKU.link
A_calving(end,end) = A_calving(end-1,end)-4.320; 

M_calving = [M_calving;M_calving(end,:)];
M_calving(end,132) = M_calving(end-1,132)-4.320*ronne_mass_per_area; 
M_calving(end,end) = M_calving(end-1,end)-4.320*ronne_mass_per_area; 

M_calving_err = [M_calving_err;M_calving_err(end,:)];

%%

[xc,yc] = polycenter(D.x,D.y);
[latc,lonc] = ps2ll(xc,yc);
for k = 1:181
   tmpx = D.x{k}; 
   tmpy = D.y{k}; 
   P(k) = polyshape(tmpx,tmpy); 
end
A = area(P); 

col = mat2rgb(lonc,cmocean('phase'),[-1 1]*180); 

%%
HMdiff = dHdt*(year(end)-year(1)); 
CMdiff = (M_calving(end,:)- M_calving(1,:))'; 
CMerr =abs(M_calving_err(end,:)- M_calving_err(1,:))';

HMdiff(end-1:end) = []; 
CMdiff(end-1:end) = []; 

% 
% CMdiff = sqrt(abs(CMdiff)).*sign(CMdiff); 
% HMdiff = sqrt(abs(HMdiff)).*sign(Mdiff); 


%figure
%hdot = plot(HMdiff,CMdiff,'.','color',rgb('black'));
Asc = A.^.5; 



%Flannery Appearance Compensation case:
% Pmin=2.5; % minimum point size 
% Pj = 1.0083 * (Asc./min(Asc)).^0.5716 * Pmin;

[~,ind] = sort(Asc,'descend'); 

figure('pos',[10 400 353*1.6 265*1.6])
hold on
for k = 1:length(ind)
   hdots(k) = scatter(HMdiff(ind(k)),CMdiff(ind(k)),90*Asc(ind(k))/max(Asc),col(ind(k),:),'filled');
   %hdots(k) = scatter(HMdiff(ind(k)),CMdiff(ind(k)),Pj(ind(k)),col(ind(k),:),'filled');
end
% for k = 1:length(ind)
%    hdots(k) = scatter(HMdiff(ind(k)),CMdiff(ind(k)),10,col(ind(k),:),'filled');
% end
%    
%hdots = scatter(HMdiff,CMdiff,300*(A.^(.75))/(max(A.^(.75))),lonc,'filled');
%hdots = scatter(HMdiff,CMdiff,30,lonc,'filled');
hold on
% for k = 1:181
%    hdot(k) = plot(HMdiff(k) +dHdt_err(k)*(year(end)-year(1))*[-1 1],CMdiff(k)*[1 1],'color',col(k,:)); 
% end
% for k = 1:181
%    hdot2(k) = plot(HMdiff(k)*[1 1],CMdiff(k) + CMerr(k)*[-1 1],'color',col(k,:)); 
% end

% for k = 1:181
%    hdot(k) = plot(HMdiff(k) +dHdt_err(k)*(year(end)-year(1))*[-1 1],CMdiff(k)*[1 1],'color','k'); 
% end
% for k = 1:181
%    hdot2(k) = plot(HMdiff(k)*[1 1],CMdiff(k) + CMerr(k)*[-1 1],'color','k'); 
% end

caxis([-1 1]*180)
cmocean phase
sidelabel = 2040;
txt(1) = text(-sidelabel,0,{'T';'H';'I';'N';'N';'E';'D'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(2) = text(sidelabel,0,{'T';'H';'I';'C';'K';'E';'N';'E';'D'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(3) = text(0,-sidelabel,{'RETREATED'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');
txt(4) = text(0,sidelabel,{'ADVANCED'},'color','w',...
   'fontweight','bold','fontsize',12,'horiz','center');

try
   uistack(hdot(:),'top'); 
uistack(hdot2(:),'top'); 
catch
uistack(flipud(hdots(:)),'top'); 
end
ind = hypot(HMdiff,CMdiff)>110;
tmp = HMdiff; 
ind(124)=true; 
ind(83) = true; 
ind(155) = true; 
tmp(~ind) = nan; 
lt = text(tmp,CMdiff,D.name,'horiz','center',...
   'vert','bot','fontsize',5,'color',rgb('black'));

set(lt([84:90 121 135 170 124 83]),'vert','top')
set(lt(22),'horiz','left','vert','mid','string',' Brunt Stancomb ') % brunt
set(lt(168),'horiz','right','vert','mid','string',' West ') % west
set(lt(49),'horiz','left','vert','bot','string',' Fimbul ') % fimbul
set(lt(105),'horiz','center','vert','top','string',' Mertz ') % mertz
set(lt(1),'horiz','center','vert','top') % abbot
set(lt(60),'horiz','right','vert','top') % george
set(lt(39),'horiz','right','vert','mid','string',' Dotson ') % 
set(lt(113),'horiz','right','vert','top','string',' Ninnis') % 
set(lt(84),'horiz','right','vert','mid','string',' LarsenA ') % 
set(lt(134),'horiz','left','vert','top','string',' Ross East ') % 
set(lt(155),'horiz','left','vert','mid','string',' Totten ') % 
set(lt(124),'horiz','left','vert','top','string',' Prince Harald ') % fimbul
set(lt(146),'string','') %stange was getting in the way
axis equal
%axis([-1 1 -1 1]*1.01*max(abs(axis)))

axis((sidelabel+100)*[-1 1 -1 1])

% xlabel('Mass change due to thinning or thickening (Gt)','fontsize',7)
% ylabel('Mass change due to retreat or advance (Gt)','fontsize',7)
title 'Ice shelf mass change since 1997 (Gt)' 

hold on
p(1)=hline(0,'color',.83*[1 1 1]);
p(2)=vline(0,'color',.81*[1 1 1]);
%p(3) = plot(xlim,xlim,'color',rgb('light gray')); 
uistack(p,'bottom')
set(gca,'fontsize',7,'yticklabelrotation',90)

ax = axis; 
box on
%set(gca,'color',rgb('pastel pink')); 
cm = crameri('-vik',1001); % broc, cork, vik
set(gca,'color',cm(501-60,:)); 

hp = patch([ax(1) ax(2) 0 ax(1)],[ax(3) ax(3) 0 ax(3)],'b'); 
%hp.FaceColor = rgb('pastel blue'); 
hp.FaceColor = cm(501+60,:); 

hp.EdgeColor = 'none'; 
uistack(hp,'bottom') 

hp(2) = patch([ax(1) 0 ax(2) ax(1)],[ax(4) 0 ax(4) ax(4)],'b'); 
%hp(2).FaceColor = rgb('pastel blue'); 
hp(2).FaceColor = cm(501+60,:) ;

hp(2).EdgeColor = 'none'; 
uistack(hp(2),'bottom') 

% export_fig calving_vs_thinning_masschange.png -r600 -p0.01 

%

gp = plotboxpos(gca);
axes('position',[gp(1)+(.75)*gp(3) gp(2)+(.75)*gp(4) .25*gp(3) .25*gp(4)])
 
hold on
for k = 1:181
plot(P(k),'facecolor',col(k,:),'facealpha',1,'edgecolor','none')
end
% bedmachine('gl','linewidth',0.2,'color',0.5*[1 1 1])
% axis image
%hant = antbounds('coast','polyshape','facecolor','w','facealpha',1,'edgecolor',.5*[1 1 1],'linewidth',0.2);
hant = antbounds('coast','polyshape','facecolor','w','facealpha',1,'edgecolor','none');
uistack(hant,'bottom'); 
axis tight off
bedmachine('coast','color',0.3*[1 1 1],'linewidth',0.1)
bedmachine('gl','color',0.3*[1 1 1],'linewidth',0.1)
gp2 = plotboxpos(gca); 
set(gca,'xcolor','none','ycolor','none','pos',[gp(1)+gp(3)-gp2(3) gp(2)+gp(4)-gp2(4) gp2(3) gp2(4)])
 
%exportgraphics(gcf,'calving_vs_thinning_masschange.pdf'); 
%export_fig calving_vs_thinning_masschange.pdf -painters 

%%

A1 = A_calving(1,1:181); 
dA = A_calving(25,1:181)-A_calving(1,1:181); 
[A1,ind] = sort(A1); 
dA = dA(ind); 

for k=1:181
   N_losers_percent(k) = 100*sum(dA(1:k)<0)/k;
   A_loss_below(k) = mean(100*dA(1:k)./A1(1:k)); 
end




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

