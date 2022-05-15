
% This script makes a simple cartoon to show how we use 2D interpolation 
% to get a calving flux estimate of a gridded dataset. 
% Chad A. Greene, NASA JPL, October 2021. 

%% Preferences: 

gridcol = 0.5*[1 1 1];     % color of the grid
newcol = hex2rgb('c76084');% color of the vectors and new grid box
newcol2 = hex2rgb('959b44'); % color of just the part of the new grid cell that has calved. 

%% Create the grid: 

N = 5; % to create a NxN grid. You pretty much don't want to chane this value, unless you also want to adjust the indices of x0,y0 and x1,y1
x = 1:N; 
y = 1:N; 

% Create the ice mask: 
ice = flipud(tril(true(N))); 
ice(4:5,:) = false; 

% Prescribe displacements: 
[X,Y] = meshgrid(x,y); 
dX = 0.2*ones(size(X)); 
dY = 0.1*ones(size(X)); 

dX(~ice) = NaN; 
dY(~ice) = NaN; 

% Center coordinates of the center pixel (the star of the show!):
x0 = x(3); 
y0 = y(3); 
x1 = x0+dX(3,3); 
y1 = y0+dY(3,3); 

%% Make a figure

figure('pos',[100 100 330 237])
imagescn(x,y,ice)
cmocean ice
caxis([-0.1 1]) % lightens up the ocean a tiny bit
hold on
plot(X(:),Y(:),'.','color',gridcol) % dots at grid cell centers. 

vline(0.5:(N+1),':','color',gridcol)
hline(0.5:(N+1),':','color',gridcol)
axis off

hdot = plot(x1,y1,'.','color',newcol); 

% Outline the new location of the center pixel: 
hbox = patch(x1+0.5*[-1 1 1 -1 -1],y1+0.5*[-1 -1 1 1 -1],'r');
hbox.FaceColor = newcol; 
hbox.FaceAlpha = 0.05; 
hbox.EdgeColor = newcol; 
hbox.LineStyle = ':';

% Outline just the calved portion of the center pixel: 
hbox2 = patch([x1-.5 x0+.5 x0+.5 x1+.5 x1+.5 x1-.5 x1-.5],...
   [y0+.5 y0+.5 y1-.5 y1-.5 y1+.5 y1+.5 y0+.5],'r');
hbox2.FaceColor = newcol2; 
hbox2.FaceAlpha = 0.4; 
hbox2.EdgeColor = newcol2; 
hbox2.LineStyle = '-';

% Show displacements: 
q = quiver(X,Y,dX,dY,'r'); 
q.Color = newcol; 
q.AutoScale = 'off'; 

% Label things: 
txt(1) = text(2,3.5,' ocean','color',gridcol,'vert','bottom','horiz','center','fontsize',7);
txt(2) = text(2,3.5,' ice','color',gridcol,'vert','top','horiz','center','fontsize',7);
txt(3) = text(x0,y0,'({\itx}_0,{\ity}_0)','color',gridcol,'vert','top','horiz','center','fontsize',6,'fontname','times new roman');
txt(4) = text(x1,y1,'({\itx}_1,{\ity}_1)','color',newcol,'vert','bot','horiz','center','fontsize',6,'fontname','times new roman');
txt(5) = text(x1,y0+0.5+dY(3,3)/2,'calved ice','color',newcol2,'vert','middle','horiz','center','fontsize',7);

% Zoom in and set the aspect ratio: 
axis([1 5 1 5]+[.3 -.3 .6 -.6])
daspect([1 1 1])

% export_fig calving_flux_calculation_cartoon.png -r600
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/fig_ED5.eps','contenttype','vector')
