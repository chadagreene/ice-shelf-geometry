

%% Load data 

load('calving_flux_timeseries.mat') 
D = load('iceshelves_2008_v2.mat');

%%

figure
subplot(2,1,1) 
plot(year,(A_calving-A_calving(1,:))/1e9)
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


figure
plot(year,A_calving(:,132))
hold on
polyfit(year(5:end),A_calving(5:end,132),1)
ft = polyfit(year(5:end),A_calving(5:end,132),1)
plot(2005:2100,polyval(ft,2005:2100))
yline(A_calving(1,132))

figure
plot(year,M_calving(:,132))
hold on
ftm = polyfit(year(5:end),M_calving(5:end,132),1)
plot(2005:2100,polyval(ftm,2005:2100))
yline(M_calving(1,132))
title mass


