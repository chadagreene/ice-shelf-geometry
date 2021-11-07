

%% Load data 

load('calving_flux_timeseries.mat') 
D = load('iceshelves_2008_v2.mat');

% Convert to thousands of square km: 
A_calving = A_calving*1e-9; 

%% Figure out the mass of Ronne's 2021 calving event
% 
ronne_mass_per_area = diff(A_calving(:,132))\diff(M_calving(:,132));

figure('pos',[101.00        636.00        292.00        224.00]) 

hold on
set(gca,'fontsize',7) 
plot(diff(A_calving(:,132)),diff(M_calving(:,132)),'.','markersize',8)
hold on
axis([-11 1 -2550 200])
polyplot(diff(A_calving(:,132)),diff(M_calving(:,132)),1)
box off
xlabel 'observed change in area (thousands of km^2)'
ylabel 'observed change in mass (Gt)' 
ntitle('Ronne Ice Shelf','fontsize',7)
round(ronne_mass_per_area)

pv = polyfit(diff(A_calving(:,132)),diff(M_calving(:,132)),1);
text(-8,polyval(pv,-8),[' ',num2str(round(ronne_mass_per_area)),' Gt per thousand km^2'],'fontsize',6)

plot(-4.32,polyval(pv,-4.32),'k+')
text(-4.32,polyval(pv,-4.32),' Iceberg A76 ','vert','top','horiz','left','fontsize',6); 

% export_fig ronne_area_to_mass.png -pdf -r600 -painters

%% Add the 2021 Ronne calving to the time series: 

year = [year 2021.4]; 
A_calving = [A_calving;NaN(1,183)];
A_calving(end,132) = A_calving(end-1,132)-4.320; % Iceberg A-76 https://www.esa.int/ESA_Multimedia/Images/2021/05/Meet_the_world_s_largest_iceberg#.YKTfep3wFKU.link
A_calving(end,end) = A_calving(end-1,end)-4.320; 

M_calving = [M_calving;NaN(1,183)];
M_calving(end,132) = M_calving(end-1,132)-4.320*ronne_mass_per_area; 
M_calving(end,end) = M_calving(end-1,end)-4.320*ronne_mass_per_area; 

M_calving_err = [M_calving_err;M_calving_err(end,:)];

yearc = year(1:end-1)+diff(year)/2; 

%%

figure
subplot(2,1,1) 
plot(year,(A_calving-A_calving(1,:)))
box off
axis tight
ax(1) = gca; 
ylabel 'area change (thousands of km^2') 

subplot(2,1,2)
plot(year,M_calving-M_calving(1,:))
box off
axis tight
ax(2) = gca; 
ylabel 'mass change (Gt') 
linkaxes(ax,'x')


%% Time series of past and future calving

col = hex2rgb({'#50b47b';
'#b94b75';
'#7f63b8';
'#98a441';
'#ba6437'});

extrapyr = 2021.5:2040; 
ind = 8;
figure('pos',[19.00        553.00        639.00        271.00]) 
hold on
set(gca,'fontsize',7) 
fs1 = 5; 
offset = 0; 
fw = 'bold'; 

% Ronne A-38 October 1998, and May 2000 A-43; also May 2021 A-76 4320 km^2 
ft = polyfit(year(11:end-1),A_calving(11:end-1,132),1);
tmp = A_calving(:,132) - (A_calving(1,132) +ft(1)*1);
yr3 = [year(1) 1998.8 year(2:end)]; 
tmp3 = [tmp(1) 0 tmp(2:end)']; 
plot(yr3,tmp3,'.-','color',col(1,:),'linew',1,'markersize',10)  
hold on
plot(extrapyr,polyval(ft,extrapyr) - (A_calving(1,132) +ft(1)*1)-4.32,':','color',col(1,:),'linew',1)
text(year(ind),tmp(ind)+offset,'Ronne','horiz','left','vert','top','fontsize',6,'fontweight',fw,'color',col(1,:))
%text(1998.8,0,{'A38-39: 8078 km^2';'Oct 13, 1998';'\it{Lazzara, 1999}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(1,:))
text(mean(yr3(2:3)),mean(tmp3(2:3)),{'A38-39 '},'vert','mid','horiz','right','fontsize',fs1,'color',col(1,:))
text(mean(yr3(end-1:end)),mean(tmp3(end-1:end)),{' A76 '},'vert','mid','horiz','left','fontsize',fs1,'color',col(1,:))
text(mean(yr3(3:4)),mean(tmp3(3:4)),{' A43-44  '},'vert','mid','horiz','right','fontsize',fs1,'color',col(1,:))


% Ross East: C-19 May 2002 
ft = polyfit(year(6:end-1),A_calving(6:end-1,134),1); 
tmp = A_calving(:,134) - (A_calving(3,134) +ft(1)*2/12);
yr3 = [year(1:5)  year(6:end)]; 
tmp3 = [tmp(1:5)' tmp(6:end)']; 
plot(yr3,tmp3,'.-','color',col(2,:),'linew',1,'markersize',10)
plot(extrapyr,polyval(ft,extrapyr) - (A_calving(3,134) +ft(1)*2/12),':','color',col(2,:),'linew',1)
text(year(ind),tmp(ind)+offset,'Ross East','horiz','left','vert','top','fontsize',6,'fontweight',fw,'color',col(2,:))
%text(2002.36+.2,0,{'C19: 6368 km^2';'May 2002';'\it{Budge & Long, 2018}'},'vert','top','horiz','left','fontsize',fs1,'color',col(2,:))
text(mean(yr3(5:6)),mean(tmp3(5:6)),{' C18-19'},'vert','top','horiz','LEFT','fontsize',fs1,'color',col(2,:))
text(mean(yr3(3:4)),mean(tmp3(3:4)),{' C16'},'vert','mid','horiz','LEFT','fontsize',fs1,'color',col(2,:))

% Western Ross: B-15 March 2000 
ft = polyfit(year(4:end-1),A_calving(4:end-1,135),1);
tmp = A_calving(:,135)-(A_calving(1,135)+ft(1)*2.5);
yr3 = [year(1) 2000.18 year(2:end)]; 
tmp3 = [tmp(1) 0 tmp(2:end)']; 
plot(yr3,tmp3,'.-','color',col(3,:),'linew',1,'markersize',10)
plot(extrapyr,polyval(ft,extrapyr)-(A_calving(1,135)+ft(1)*2.5),':','color',col(3,:),'linew',1)
text(year(ind),tmp(ind)+offset,'Ross West','horiz','left','vert','top','fontsize',6,'fontweight',fw,'color',col(3,:))
% B-9 October 1987, 4540 km^2 keys1990calving
plot([1987.7 1987.8 year(1)],[0 -4.54 tmp(1)],':.','color',col(3,:),'linew',1,'markersize',10)
%text(1987.7,0,{'B9: 4540 km^2';'Oct 1987';'\it{Keys et al., 1990}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(3,:))
%text(2000.25,0,{'B15-18: 13,135 km^2';'Mar-Apr 2000';'\it{Lazzara et al., 1999}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(3,:))
text(1987.75,-4.54/2,{' B8-9'},'vert','mid','horiz','left','fontsize',fs1,'color',col(3,:))
text(mean(yr3(2:3)),mean(tmp3(2:3)),{'B15-18 '},'vert','bot','horiz','right','fontsize',fs1,'color',col(3,:))

% Filchner 1986 calving event created A22, A23, and maybe A24 but I'm not sure. 
ft = polyfit(year(4:end-1),A_calving(4:end-1,48),1);
FilchnerGrowthRate = ft(1)*1000
tmp = A_calving(:,48)-(A_calving(1,48)-ft(1)*11)-11.5; % 11,500 km^2 calving event in 1986 according to Ferringo1987substantial 
plot(year,tmp,'.-','color',col(4,:),'linew',1,'markersize',10)
ft = polyfit(year(11:end-1),tmp(11:end-1),1);
plot(extrapyr,polyval(ft,extrapyr),':','color',col(4,:),'linew',1)
text(year(ind),tmp(ind)+offset,'Filchner','horiz','left','vert','top','fontsize',6,'fontweight',fw,'color',col(4,:))
plot([1986 1986.1 year(1)],[0 -11.5 tmp(1)],':.','color',col(4,:),'linew',1,'markersize',10)
%text(1986,0,{'11,500 km^2';'Apr to Nov 1986';'\it{Ferrigno & Gould, 1987}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(4,:))
text(1986.05,-11.5/2,{'Icebergs ';'A22-24 '},'vert','bot','horiz','right','fontsize',fs1,'color',col(4,:))

% Amery 10,000 km^2 in 1963 or 1964 (Fricker2002iceberg)
ft = polyfit(year(1:end-3),A_calving(1:end-3,10),1);
AmeryGrowthRate = ft(1)*1000
tmp = A_calving(:,10)-(A_calving(1,10)-ft(1)*34.25)-10; % 10 thousand km^2 event 34.25 years before 1997 
%plot(year,tmp,'.-','color',col(5,:),'linew',1,'markersize',10)

yr3 = [year(1:22) 2019.73 year(23:end)]; 
tmp3 = [tmp(1:22)' tmp(22)+ft(1)*.5 tmp(23:end)'];
tmp4 = tmp3; 
plot(yr3,tmp3,'.-','color',col(5,:),'linew',1,'markersize',10)
tmp3 = polyval(ft,extrapyr); 
plot(extrapyr,tmp3-tmp3(1)+tmp(end-1),':','color',col(5,:),'linew',1)
text(year(ind),tmp(ind)+offset,'Amery','horiz','left','vert','top','fontsize',6,'fontweight',fw,'color',col(5,:))
plot([1963.5 1963.55 year(1)],[0 -10 tmp(1)],':.','color',col(5,:),'linew',1,'markersize',10)
%text(1963.5,0,{'10,000 km^2 calved in';'late 1963 or early 1964';'\it{Fricker et al., 2002}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(5,:))
%text(2019.75,0,{'D28: 1636 km^2';'Sept 25, 2019';'\it{Walker et al., 2021}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(5,:))
text(1963.525,-5,{' 10,000 km^2 calving event'},'vert','bot','horiz','left','fontsize',fs1,'color',col(5,:))
text(1963.525,-5,{' (\it{Fricker et al., 2002})'},'vert','top','horiz','left','fontsize',fs1,'color',col(5,:))
text(mean([2019.73 year(23)]),mean(tmp4(23:24)),{' D28'},'vert','mid','horiz','left','fontsize',fs1,'color',col(5,:))

box off 
axis tight
hyline = plot([1962 max(extrapyr)],[0 0],'color',rgb('gray'),'linewidth',0.4); 
uistack(hyline,'bottom') 
ylabel 'area w.r.t. maximum known extents (thousands of km^2)'

axis([1962 max(extrapyr) -18 1.9])
set(gca,'xtick',1965:5:2040)
% export_fig calving_past_and_future.png -r600 -painters
% export_fig calving_past_and_future.pdf -r600 -painters
%%

fs1 = 5; 


rng('default')
rng(2) % seed the random number generator for repeatable colors 
% col = rand(181,3); 

col = double(lab2rgb([80*rand(181,1)+5 128*2*(rand(181,1)-0.5) 128*2*(rand(181,1)-0.5)],'OutputType','uint8'))/255;

A = A_calving-A_calving(1,:);

figure('pos',[12.00        502.00        324.00        400.00]) 
hold on
set(gca,'fontsize',7)
h = plot(year,A(:,1:181)); 
for k = 1:181
   h(k).Color = col(k,:); 
   if k==132   
      txt(k)=text(year(end),A(end,k),[' ',D.name{k}],'fontsize',fs1,'color',col(k,:),'vert','middle');
   else
      if (abs(A(end-1,k))>A(end-1,140) | k==1) % if the magnitude is greater than shackleton (because shackleton is where labels start to overlap)
         txt(k)=text(year(end-1),A(end-1,k),[' ',D.name{k}],'fontsize',fs1,'color',col(k,:),'vert','middle');
      end
   end
end

big5 = false(1,183); 
big5([10 48 132 134 135]) = true; 
notbig5 = ~big5;
notbig5(end) = false; 
%plot(year,sum(A(:,big5),2,'omitnan'))
plot(year,A(:,end)-sum(A(:,big5),2,'omitnan'))
plot(year,A(:,end)-sum(A(:,big5),2,'omitnan'))
plot(year,sum(A(:,notbig5),2,'omitnan')); 

% Fix a couple of labels: 
%txt(10).VerticalAlignment = 'top'; % amery
%txt(1).VerticalAlignment = 'bottom'; % abbot
txt(22).VerticalAlignment = 'bottom'; % brunt
txt(31).VerticalAlignment = 'bottom'; % cook
txt(86).VerticalAlignment = 'top'; % larsen C
%txt(168).VerticalAlignment = 'top'; % west

hold on
plot(year,A(:,end),'k-','linewidth',1) 
text(year(end),A(end,end),' Antarctica','fontsize',fs1,'fontweight','bold')
axis tight
box off
xlim([year(1) 2024]); 
pl = plot(xlim,[0 0],'color',rgb('gray'),'linewidth',0.4); 
uistack(pl,'bottom')
ylabel 'change in ice shelf area (thousands of km^2)'

set(gca,'outerpos',[0 0 .9 1])

% export_fig calving_timeseries.png -r600 -painters
% export_fig calving_timeseries.pdf -r600 -painters

%% Area chart 

fs1 = 5; 


rng('default')
rng(2) % seed the random number generator for repeatable colors 
% col = rand(181,3); 

col = double(lab2rgb([80*rand(181,1)+5 128*2*(rand(181,1)-0.5) 128*2*(rand(181,1)-0.5)],'OutputType','uint8'))/255;

A = A_calving-A_calving(1,:);

figure('pos',[12.00        502.00        307.00        379.00]) 
hold on
set(gca,'fontsize',7)

Atmp = A(:,1:182); 
D.name{182} = '';
[Atmp,ind] = sort(sum(Atmp,'omitnan')); 
%[Atmp,ind] = sort(min(Atmp(end-1:end,:))); 
As = A(:,ind); 
isn = isnan(As(end,:)); 
As(end,isn) = As(end-1,isn); 

A_lo = fliplr(As(:,Atmp<0)); 
A_hi = As(:,Atmp>0); 
h_lo = area(year,A_lo,'edgecolor','none');
h_hi = area(year,A_hi,'edgecolor','none');

A_loc = [zeros(25,1) cumsum(A_lo,2)];
A_hic = [ zeros(25,1) cumsum(A_hi,2)];

name_2 = D.name(ind); 
name_lo = name_2(Atmp<0); 
name_lo = flipud(name_lo); 
name_hi = name_2(Atmp>0); 


cm = crameri('acton',22); 
cmind = 3+[1:5:19 2:5:19 3:5:19 4:5:19 5:5:19];
cm = repmat(cm(cmind,:),10,1); 

for k = 1:length(h_lo)
   h_lo(k).FaceColor = cm(length(h_lo)-k+1,:); 
end
plot(year,A_loc(:,end),'color',cm(1,:),'linewidth',1); 

cm = crameri('oslo',25); 
cm(1:2,:) = []; 
cmind = 2+[1:5:19 2:5:19 3:5:19 4:5:19 5:5:19];
cm = repmat(cm(cmind,:),10,1); 
for k = 1:length(h_hi)
   h_hi(k).FaceColor = cm(k,:); 
end


A_loc = [zeros(25,1) cumsum(A_lo,2)]; 
for k = 91:length(h_lo)
   try
   [xc,yc] = polycenter([year fliplr(year)]',[A_loc(:,k);flipud(A_loc(:,k+1))]); 
   text(xc,yc,name_lo{k},'horiz','center','vert','mid','fontsize',5,'color',hex2rgb('ffffff'));
   end
end

A_hic = [ zeros(25,1) cumsum(A_hi,2)]; 
for k = 67:length(h_hi)
   try
   [xc,yc] = polycenter([year fliplr(year)]',[A_hic(:,k);flipud(A_hic(:,k+1))]); 
   text(xc,yc,name_hi{k},'horiz','center','vert','mid','fontsize',5,'color',hex2rgb('ffffff'));
   end
end
box off
axis tight
set(gca,'yaxislocation','left')


colorcase = 63; 
switch colorcase 
   case 0
      linecol = [0 0 0]; % color of the Antarctica line
      shadowcol = .8*[1 1 1]; 
   case 3
      linecol = [1 1 1]; % color of the Antarctica line
      shadowcol = 0*[1 1 1]; 

   case 5
      shadowcol = [0 0 0];
      linecol = hex2rgb('88d895'); 
   case 6
      shadowcol = [0 0 0];
      linecol = hex2rgb('daba63'); 
   case 62
      shadowcol = [0 0 0];
      linecol = hex2rgb('f7e78f'); 
   case 63
      shadowcol = [0 0 0];
      linecol = hex2rgb('f0d079'); 
   case 7
      linecol = rgb('gold');
      shadowcol = [0 0 0]; 
   case 8
      shadowcol = [0 0 0];
      linecol = hex2rgb('bbe19d'); 

      
      
   otherwise
      error 'unknown colorcase'
end


plot(year,A(:,end),'color',shadowcol,'linewidth',2)
plot(year(end),A(end,end),'.','linewidth',2,'color',shadowcol,'markersize',8)
plot(year(1),0,'.','linewidth',2,'color',shadowcol,'markersize',8)
plot(year,A(:,end),'color',linecol,'linewidth',1.5)
plot(year(end),A(end,end),'.','linewidth',1.5,'color',linecol,'markersize',7)
plot(year(1),0,'.','linewidth',1.5,'color',linecol,'markersize',7)
%text(year(21),A(21,end),'Antarctica','fontsize',6,'fontweight','bold','horiz','center','vert','top');
%text(year(17),A(17,end),'Antarctica','fontsize',6,'fontweight','bold','horiz','center','vert','top');
%text(year(16),A(16,end),'Antarctica','fontsize',6,'fontweight','bold','horiz','center','vert','bot');
%text(year(end),A(end,end),'Antarctica ','fontsize',6,'fontweight','bold','horiz','right','vert','mid','color',linecol);

% text(year(end),A(end,end),'Antarctica ','fontname','Helvetica Neue','fontsize',7,'fontweight','bold','horiz','right','vert','mid','color',linecol);
% text(year(end),A(end,end),'Antarctica ','fontname','HelveticaNeue LT 75 BdOutline','fontweight','bold','fontsize',7,'horiz','right','vert','mid','color',shadowcol);
drawnow
textborder(year(end),A(end,end),'Antarctica ',linecol,shadowcol,'fontname','Helvetica','fontsize',6.5,'fontweight','bold','horiz','right','vert','mid')

axcol = 0.5*[1 1 1]; 
xl = xlim; 
hl = plot(xl,repmat((-50:10:20)',1,2),'-','color',axcol); 
uistack(hl,'bottom')
txt = text(xl(1)*ones(length(-50:10:20),1),(-50:10:20)',num2str((-50:10:20)'),'horiz','left','vert','top','fontsize',7,'color',axcol);
uistack(txt(:),'bottom')
txt(end).String = '+20\times1000 km^2'; 
txt(end-1).String = '+10'; 
text(xl(1),20,'Cumulative area change since 1997','fontsize',7,'horiz','left','vert','bot','fontweight','bold','color',.1*[1 1 1])
set(gca,'ycolor','none','xcolor',axcol)

% export_fig cumulative_area_change.png -r600 -painters -p0.01
%exportgraphics(gca,'cumulative_area_change.pdf')
% export_fig test.png -r600 -painters -p0.01
%% (OLD) Time series of past and future calving 
% 
% col = hex2rgb({'#50b47b';
% '#b94b75';
% '#7f63b8';
% '#98a441';
% '#ba6437'});
% 
% extrapyr = 2021.5:2040; 
% ind = 10;
% %figure('pos',[19.00        553.00        560.00        337.00]) 
% figure('pos',[19.00        553.00        794.00        337.00]) 
% 
% % Ronne A-38 October 1998, and May 2000 A-43; also May 2021 A-76 4320 km^2 
% ft = polyfit(year(11:end-1),A_calving(11:end-1,132),1);
% tmp = A_calving(:,132) - (A_calving(1,132) +ft(1)*1);
% yr3 = [year(1) 1998.8 year(2:end)]; 
% tmp3 = [tmp(1) 0 tmp(2:end)']; 
% plot(yr3,tmp3,'.-','color',col(1,:),'linew',1,'markersize',10)  
% hold on
% plot(extrapyr,polyval(ft,extrapyr) - (A_calving(1,132) +ft(1)*1)-4.32,':','color',col(1,:),'linew',1)
% text(year(ind),tmp(ind),'Ronne','horiz','left','vert','top','fontsize',8,'color',col(1,:))
% %text(1998.8,0,{'A38-39: 8078 km^2';'Oct 13, 1998';'\it{Lazzara, 1999}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(1,:))
% text(1998.8,0,{'A38-39'},'vert','bot','horiz','right','fontsize',fs1,'color',col(1,:))
% text(yr3(end-1),tmp3(end-1),{' A76 '},'vert','top','horiz','left','fontsize',fs1,'color',col(1,:))
% text(yr3(3),tmp3(3),{' A43-44 '},'vert','top','horiz','right','fontsize',fs1,'color',col(1,:))
% 
% 
% % Ross East: C-19 May 2002 
% ft = polyfit(year(6:end-1),A_calving(6:end-1,134),1); 
% tmp = A_calving(:,134) - (A_calving(3,134) +ft(1)*2/12);
% yr3 = [year(1:5)  year(6:end)]; 
% tmp3 = [tmp(1:5)' tmp(6:end)']; 
% plot(yr3,tmp3,'.-','color',col(2,:),'linew',1,'markersize',10)
% plot(extrapyr,polyval(ft,extrapyr) - (A_calving(3,134) +ft(1)*2/12),':','color',col(2,:),'linew',1)
% text(year(ind),tmp(ind),'Ross East','horiz','left','vert','top','fontsize',8,'color',col(2,:))
% %text(2002.36+.2,0,{'C19: 6368 km^2';'May 2002';'\it{Budge & Long, 2018}'},'vert','top','horiz','left','fontsize',fs1,'color',col(2,:))
% text(yr3(5),tmp3(5),{' C18-19'},'vert','mid','horiz','LEFT','fontsize',fs1,'color',col(2,:))
% text(yr3(3),tmp3(3),{' C16'},'vert','top','horiz','LEFT','fontsize',fs1,'color',col(2,:))
% 
% 
% % Western Ross: B-15 March 2000 
% ft = polyfit(year(4:end-1),A_calving(4:end-1,135),1);
% tmp = A_calving(:,135)-(A_calving(1,135)+ft(1)*2.5);
% yr3 = [year(1) 2000.18 year(2:end)]; 
% tmp3 = [tmp(1) 0 tmp(2:end)']; 
% plot(yr3,tmp3,'.-','color',col(3,:),'linew',1,'markersize',10)
% plot(extrapyr,polyval(ft,extrapyr)-(A_calving(1,135)+ft(1)*2.5),':','color',col(3,:),'linew',1)
% text(year(ind),tmp(ind),'Ross West','horiz','left','vert','top','fontsize',8,'color',col(3,:))
% % B-9 October 1987, 4540 km^2 keys1990calving
% plot([1987.7 1987.8 year(1)],[0 -4.54 tmp(1)],':.','color',col(3,:),'linew',1,'markersize',10)
% %text(1987.7,0,{'B9: 4540 km^2';'Oct 1987';'\it{Keys et al., 1990}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(3,:))
% %text(2000.25,0,{'B15-18: 13,135 km^2';'Mar-Apr 2000';'\it{Lazzara et al., 1999}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(3,:))
% text(1987.7,0,{'B8-9'},'vert','bot','horiz','center','fontsize',fs1,'color',col(3,:))
% text(2000.25,0,{'B15-18'},'vert','bot','horiz','center','fontsize',fs1,'color',col(3,:))
% 
% % Filchner 1986 calving event created A22, A23, and maybe A24 but I'm not sure. 
% ft = polyfit(year(4:end-1),A_calving(4:end-1,48),1);
% tmp = A_calving(:,48)-(A_calving(1,48)-ft(1)*11)-11.5; % 11,500 km^2 calving event in 1986 according to Ferringo1987substantial 
% plot(year,tmp,'.-','color',col(4,:),'linew',1,'markersize',10)
% ft = polyfit(year(11:end-1),tmp(11:end-1),1);
% plot(extrapyr,polyval(ft,extrapyr),':','color',col(4,:),'linew',1)
% text(year(ind),tmp(ind),'Filchner','horiz','left','vert','top','fontsize',8,'color',col(4,:))
% plot([1986 1986.1 year(1)],[0 -11.5 tmp(1)],':.','color',col(4,:),'linew',1,'markersize',10)
% %text(1986,0,{'11,500 km^2';'Apr to Nov 1986';'\it{Ferrigno & Gould, 1987}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(4,:))
% text(1986,0,{'Icebergs A22-24 '},'vert','bot','horiz','right','fontsize',fs1,'color',col(4,:))
% 
% % Amery 10,000 km^2 in 1963 or 1964 (Fricker2002iceberg)
% ft = polyfit(year(1:end-3),A_calving(1:end-3,10),1);
% tmp = A_calving(:,10)-(A_calving(1,10)-ft(1)*34.25)-10; % 10 thousand km^2 event 34.25 years before 1997 
% %plot(year,tmp,'.-','color',col(5,:),'linew',1,'markersize',10)
% 
% yr3 = [year(1:22) 2019.73 year(23:end)]; 
% tmp3 = [tmp(1:22)' tmp(22)+ft(1)*.5 tmp(23:end)']; 
% plot(yr3,tmp3,'.-','color',col(5,:),'linew',1,'markersize',10)
% tmp3 = polyval(ft,extrapyr); 
% plot(extrapyr,tmp3-tmp3(1)+tmp(end-1),':','color',col(5,:),'linew',1)
% text(year(ind),tmp(ind),'Amery','horiz','left','vert','top','fontsize',8,'color',col(5,:))
% plot([1963.5 1963.55 year(1)],[0 -10 tmp(1)],':.','color',col(5,:),'linew',1,'markersize',10)
% %text(1963.5,0,{'10,000 km^2 calved in';'late 1963 or early 1964';'\it{Fricker et al., 2002}'},'vert','bot','horiz','left','fontsize',fs1,'color',col(5,:))
% %text(2019.75,0,{'D28: 1636 km^2';'Sept 25, 2019';'\it{Walker et al., 2021}'},'vert','bot','horiz','right','fontsize',fs1,'color',col(5,:))
% text(1963.5,0,{' 10,000 km^2 calving event'},'vert','bot','horiz','left','fontsize',fs1,'color',col(5,:))
% text(1963.5,0,{' (\it{Fricker et al., 2002})'},'vert','top','horiz','left','fontsize',fs1,'color',col(5,:))
% text(2019.75,0,{'D28'},'vert','bot','horiz','center','fontsize',fs1,'color',col(5,:))
% 
% box off 
% axis tight
% hyline = yline(0,'color',rgb('gray')); 
% uistack(hyline,'bottom') 
% ylabel 'area wrt maximum known extents (thousands of km^2)'
% 
% axis([1962 max(extrapyr) -18 1.9])
% % export_fig Ross_FRIS_calving.png -r500 -p0.01 -painters
% % export_fig Ross_FRIS_calving_annotated.jpg -r600 -p0.01 -painters


%% Subfunctions 

function h = textborder(x, y, string, text_color, border_color, varargin)
%TEXTBORDER Display text with border.
%   TEXTBORDER(X, Y, STRING)
%   Creates text on the current figure with a one-pixel border around it.
%   The default colors are white text on a black border, which provides
%   high contrast in most situations.
%   
%   TEXTBORDER(X, Y, STRING, TEXT_COLOR, BORDER_COLOR)
%   Optional TEXT_COLOR and BORDER_COLOR specify the colors to be used.
%   
%   Optional properties for the native TEXT function (such as 'FontSize')
%   can be supplied after all the other parameters.
%   Since usually the units of the parent axes are not pixels, resizing it
%   may subtly change the border of the text out of position. Either set
%   the right size for the figure before calling TEXTBORDER, or always
%   redraw the figure after resizing it.
%   
%   Author: João F. Henriques, April 2010
% Heavily edited by Chad Greene for FOP tools, November 2015, including the
% following changes: 
%   * Now prints 8 background texts rather than the previous 4. 
%   * Now returns object handles 
%   * Can now handle multiple text inputs
%   * Now prints in data units rather than data units, converting to pixels, moving, 
%       then converting back to data units. This is b/c changing units was VERY slow
%       on R2014b for some reason.  Took about 35 seconds before, which is absurd.  



if nargin < 5, border_color = 'k'; end  %default: black border
if nargin < 4, text_color = 'w'; end  %default: white text

pos = getpixelposition(gca); 
xl = get(gca,'xlim'); 
xperpx = diff(xl)/pos(3); 
yl = get(gca,'xlim'); 
yperpx = diff(yl)/pos(4); 

% border around the text, composed of 8 text objects  
offsets = [xperpx yperpx].*[[0 -1; -1 0; 0 1; 1 0] ; 0.71*[ 1 1; -1 1; -1 -1; 1 -1]];

% Initialize counters for background and main text:  
cbg = 1; 
    
for n = 1:length(x)     
	for k = 1:8
		h.bg(cbg) = text(x(n)+offsets(k,1)/10, y(n)+offsets(k,2)/4, string, 'Color',border_color, varargin{:});
		
        % Increment background counter: 
        cbg = cbg+1; 
   end

end

% the actual text inside the border
h.t = text(x, y, string, 'Color',text_color, varargin{:});

if nargout==0
    clear h;
end

end