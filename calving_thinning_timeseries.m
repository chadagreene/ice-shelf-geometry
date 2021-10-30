

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