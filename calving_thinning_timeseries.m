

%% Load data 

load('calving_flux_timeseries.mat') 
D = load('iceshelves_2008_v2.mat');

% Convert to thousands of square km: 
A_calving = A_calving*1e-9; 

ronne_mass_per_area = diff(A_calving(:,132))\diff(M_calving(:,132));

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

%%

col = hex2rgb({'#50b47b';
'#b94b75';
'#7f63b8';
'#98a441';
'#ba6437'});

extrapyr = 2021.5:2040; 
ind = 10;
%figure('pos',[19.00        553.00        560.00        337.00]) 
figure('pos',[19.00        553.00        794.00        337.00]) 

% Ronne A-38 October 1998, and May 2000 A-43; also May 2021 A-76 4320 km^2 
ft = polyfit(year(11:end-1),A_calving(11:end-1,132),1);
tmp = A_calving(:,132) - (A_calving(1,132) +ft(1)*1);
yr3 = [year(1) 1998.8 year(2:end)]; 
tmp3 = [tmp(1) 0 tmp(2:end)']; 
plot(yr3,tmp3,'.-','color',col(1,:),'linew',1,'markersize',10)  
hold on
plot(extrapyr,polyval(ft,extrapyr) - (A_calving(1,132) +ft(1)*1)-4.32,':','color',col(1,:),'linew',1)
text(year(ind),tmp(ind),'Ronne','horiz','left','vert','top','fontsize',8,'color',col(1,:))
%text(1998.8,0,{'A38-39: 8078 km^2';'Oct 13, 1998';'\it{Lazzara, 1999}'},'vert','bot','horiz','right','fontsize',6,'color',col(1,:))
text(1998.8,0,{'A38-39'},'vert','bot','horiz','right','fontsize',6,'color',col(1,:))
text(yr3(end-1),tmp3(end-1),{' A76 '},'vert','top','horiz','left','fontsize',6,'color',col(1,:))
text(yr3(3),tmp3(3),{' A43-44 '},'vert','top','horiz','right','fontsize',6,'color',col(1,:))


% Ross East: C-19 May 2002 
ft = polyfit(year(6:end-1),A_calving(6:end-1,134),1); 
tmp = A_calving(:,134) - (A_calving(3,134) +ft(1)*2/12);
yr3 = [year(1:5)  year(6:end)]; 
tmp3 = [tmp(1:5)' tmp(6:end)']; 
plot(yr3,tmp3,'.-','color',col(2,:),'linew',1,'markersize',10)
plot(extrapyr,polyval(ft,extrapyr) - (A_calving(3,134) +ft(1)*2/12),':','color',col(2,:),'linew',1)
text(year(ind),tmp(ind),'Ross East','horiz','left','vert','top','fontsize',8,'color',col(2,:))
%text(2002.36+.2,0,{'C19: 6368 km^2';'May 2002';'\it{Budge & Long, 2018}'},'vert','top','horiz','left','fontsize',6,'color',col(2,:))
text(yr3(5),tmp3(5),{' C18-19'},'vert','mid','horiz','LEFT','fontsize',6,'color',col(2,:))
text(yr3(3),tmp3(3),{' C16'},'vert','top','horiz','LEFT','fontsize',6,'color',col(2,:))


% Western Ross: B-15 March 2000 
ft = polyfit(year(4:end-1),A_calving(4:end-1,135),1);
tmp = A_calving(:,135)-(A_calving(1,135)+ft(1)*2.5);
yr3 = [year(1) 2000.18 year(2:end)]; 
tmp3 = [tmp(1) 0 tmp(2:end)']; 
plot(yr3,tmp3,'.-','color',col(3,:),'linew',1,'markersize',10)
plot(extrapyr,polyval(ft,extrapyr)-(A_calving(1,135)+ft(1)*2.5),':','color',col(3,:),'linew',1)
text(year(ind),tmp(ind),'Ross West','horiz','left','vert','top','fontsize',8,'color',col(3,:))
% B-9 October 1987, 4540 km^2 keys1990calving
plot([1987.7 1987.8 year(1)],[0 -4.54 tmp(1)],':.','color',col(3,:),'linew',1,'markersize',10)
%text(1987.7,0,{'B9: 4540 km^2';'Oct 1987';'\it{Keys et al., 1990}'},'vert','bot','horiz','left','fontsize',6,'color',col(3,:))
%text(2000.25,0,{'B15-18: 13,135 km^2';'Mar-Apr 2000';'\it{Lazzara et al., 1999}'},'vert','bot','horiz','left','fontsize',6,'color',col(3,:))
text(1987.7,0,{'B8-9'},'vert','bot','horiz','center','fontsize',6,'color',col(3,:))
text(2000.25,0,{'B15-18'},'vert','bot','horiz','center','fontsize',6,'color',col(3,:))

% Filchner 1986 calving event created A22, A23, and maybe A24 but I'm not sure. 
ft = polyfit(year(4:end-1),A_calving(4:end-1,48),1);
tmp = A_calving(:,48)-(A_calving(1,48)-ft(1)*11)-11.5; % 11,500 km^2 calving event in 1986 according to Ferringo1987substantial 
plot(year,tmp,'.-','color',col(4,:),'linew',1,'markersize',10)
ft = polyfit(year(11:end-1),tmp(11:end-1),1);
plot(extrapyr,polyval(ft,extrapyr),':','color',col(4,:),'linew',1)
text(year(ind),tmp(ind),'Filchner','horiz','left','vert','top','fontsize',8,'color',col(4,:))
plot([1986 1986.1 year(1)],[0 -11.5 tmp(1)],':.','color',col(4,:),'linew',1,'markersize',10)
%text(1986,0,{'11,500 km^2';'Apr to Nov 1986';'\it{Ferrigno & Gould, 1987}'},'vert','bot','horiz','right','fontsize',6,'color',col(4,:))
text(1986,0,{'Icebergs A22-24 '},'vert','bot','horiz','right','fontsize',6,'color',col(4,:))

% Amery 10,000 km^2 in 1963 or 1964 (Fricker2002iceberg)
ft = polyfit(year(1:end-3),A_calving(1:end-3,10),1);
tmp = A_calving(:,10)-(A_calving(1,10)-ft(1)*34.25)-10; % 10 thousand km^2 event 34.25 years before 1997 
%plot(year,tmp,'.-','color',col(5,:),'linew',1,'markersize',10)

yr3 = [year(1:22) 2019.73 year(23:end)]; 
tmp3 = [tmp(1:22)' tmp(22)+ft(1)*.5 tmp(23:end)']; 
plot(yr3,tmp3,'.-','color',col(5,:),'linew',1,'markersize',10)
tmp3 = polyval(ft,extrapyr); 
plot(extrapyr,tmp3-tmp3(1)+tmp(end-1),':','color',col(5,:),'linew',1)
text(year(ind),tmp(ind),'Amery','horiz','left','vert','top','fontsize',8,'color',col(5,:))
plot([1963.5 1963.55 year(1)],[0 -10 tmp(1)],':.','color',col(5,:),'linew',1,'markersize',10)
%text(1963.5,0,{'10,000 km^2 calved in';'late 1963 or early 1964';'\it{Fricker et al., 2002}'},'vert','bot','horiz','left','fontsize',6,'color',col(5,:))
%text(2019.75,0,{'D28: 1636 km^2';'Sept 25, 2019';'\it{Walker et al., 2021}'},'vert','bot','horiz','right','fontsize',6,'color',col(5,:))
text(1963.5,0,{' 10,000 km^2 calving event'},'vert','bot','horiz','left','fontsize',6,'color',col(5,:))
text(1963.5,0,{' (\it{Fricker et al., 2002})'},'vert','top','horiz','left','fontsize',6,'color',col(5,:))
text(2019.75,0,{'D28'},'vert','bot','horiz','center','fontsize',6,'color',col(5,:))

box off 
axis tight
hyline = yline(0,'color',rgb('gray')); 
uistack(hyline,'bottom') 
ylabel 'area wrt maximum known extents (thousands of km^2)'

axis([1962 max(extrapyr) -18 1.9])
% export_fig Ross_FRIS_calving.png -r500 -p0.01 -painters
% export_fig Ross_FRIS_calving_annotated.jpg -r600 -p0.01 -painters

%%

rng(1) % seed the random number generator for repeatable colors 
col = rand(181,3); 

col = double(lab2rgb([80*rand(181,1)+5 128*2*(rand(181,1)-0.5) 128*2*(rand(181,1)-0.5)],'OutputType','uint8'))/255;

A = A_calving-A_calving(1,:);

figure
h = plot(year,A(:,1:181)); 
for k = 1:181
   h(k).Color = col(k,:); 
   if k==132   
      txt(k)=text(year(end),A(end,k),[' ',D.name{k}],'fontsize',6,'color',col(k,:),'vert','middle');
   else
      txt(k)=text(year(end-1),A(end-1,k),[' ',D.name{k}],'fontsize',6,'color',col(k,:),'vert','middle');
   end
end

% Fix a couple of labels: 
%txt(10).VerticalAlignment = 'top'; % amery
txt(22).VerticalAlignment = 'bottom'; % brunt
txt(86).VerticalAlignment = 'top'; % larsen C


hold on
plot(year,A(:,end),'k-','linewidth',2) 
text(year(end),A(end,end),' Antarctica','fontsize',6,'fontweight','bold')
axis tight
box off
ylabel 'area change (1000\timeskm^2)'
