
% This script creates complete, gridded thickness and velocity datasets, and extends
% them 100 km beyond present day coastlines. The final product is on the 

save_everything = true; 
%% Blend ITS_LIVE and Measures (Rignot) velocities: 

% Load ITS_LIVE data: 
[vxi,x,y] = itslive_data('vx'); 
vyi = itslive_data('vy'); 
wi = 1./itslive_data('v_err').^2; % weights
v_source = zeros(size(vxi),'uint8'); 

[X,Y] = meshgrid(x,y); 

% Manually remove weird bits on the edges of ice shelves in ITS_LIVE data: 
bad = inpolygon_map(x,y,[-2199380.85   -2197379.45   -2194210.57   -2187539.25   -2167942.24   -2166691.37   -2164189.63 -2162855.36   -2161687.88   -2168109.03   -2189957.61],...
   [1251348.98    1247727.13    1246103.54    1235612.66    1222499.05    1217878.07    1217253.61    1219501.66    1222124.38    1242481.69    1247852.02]) | ...
   inpolygon_map(x,y,[-1542088.64   -1541509.45   -1540133.87   -1538830.69   -1539409.88   -1539844.28   -1540133.87   -1540495.87   -1542088.64],...
   [-824542.55    -829176.08    -835402.40    -835619.59    -831637.65    -829900.07    -825990.53    -824832.14    -823528.96]) | ...
   inpolygon_map(x,y,[-1206474.73   -1215567.10   -1217521.16   -1217880.07   -1217441.41   -1214769.52   -1214011.83   -1215128.43],...
   [-1171958.09   -1176584.04   -1175626.95   -1175826.34   -1176903.07   -1179854.10   -1178338.71   -1177142.34]) | ...
   inpolygon_map(x,y,[1441247.74    1437999.46    1437187.40    1438973.94    1439623.60    1436212.91    1433127.05    1433614.30    1430690.85    1429878.78    1431340.50 1430690.85    1425331.20    1424519.13    1427280.16    1437674.64    1449206.00],...
   [-2071942.89   -2074054.27   -2076003.23   -2079413.92   -2091432.53   -2112546.30   -2113520.78   -2121966.29   -2123428.01   -2125214.56   -2127488.35 -2131061.45   -2131061.45   -2134309.72   -2136745.93   -2135609.03   -2071942.89]) ;
vxi(bad) = nan; 

% Indices of missing ITS_LIVE data: 
indi = ~isfinite(vxi) | ~isfinite(vyi) | ~isfinite(wi); 
vxi(indi) = 0; 
vyi(indi) = 0; 
wi(indi) = 0; 

% Load Measures data: 
[vxm,vym] = measures_interp('velocity',X,Y); 
wm = 1./measures_interp('error',X,Y).^2; 

% Indices of missing Measures data: 
indm = ~isfinite(vxm) | ~isfinite(vym) | ~isfinite(wm); 
vxm(indm) = 0; 
vym(indm) = 0; 
wm(indm) = 0; 

% Weighted mean velocity components: 
vx = (vxi.*wi + vxm.*wm)./(wi + wm); 
vy = (vyi.*wi + vym.*wm)./(wi + wm); 

% Replace the completely unknown values with NaNs: 
vx(indi & indm) = nan; 
vy(indi & indm) = nan; 

vx = double(vx); 
vy = double(vy); 

v_source(~indi & indm) = 1; % ITS_LIVE
v_source(indi & ~indm) = 2; % Measures v2
v_source(~indi & ~indm) = 3; % Blend of ITS_LIVE and MEasures 

clear indi indm vxi vxm vyi vym wi wm bad

%% Build a referece surface: 
% Note that BedMachine's surface is not the surface you stand on--it is the
% surface that would exist if you subtracted all FAC from the surface you stand on. In
% contrast, Bedmap2's surface is the surface you stand on. 

% Start with the BedMachine surface: 
[H_bm,x_bm,y_bm] = bedmachine_data('thickness'); 
mask_bm = bedmachine_data('mask'); 
grounded = ismember(mask_bm,[1 2 4]); 
[X_bm,Y_bm] = meshgrid(x_bm,y_bm); 

% Load REMA data: 
h_rema = wgs2gl04c(X_bm,Y_bm,rema_interp(X_bm,Y_bm,'res',200));

% Load Bedmap2 data: 
h_b2 = bedmap2_interp(X_bm,Y_bm,'surface'); 
isn = ~(h_b2>0);  
h_b2(isn) = bedmap2_interp(X_bm(isn),Y_bm(isn),'surface','nearest'); 

% load bamber: 
h_bam = wgs2gl04c(X_bm,Y_bm,bamberdem_interp(X_bm,Y_bm));
h_bam(bwdist(~isfinite(h_bam))*.5<15) = NaN; % Eliminates 15 km perimeter b/c drooping ice sheet edges in Bamber and RAMP dems

% RAMP2:
[Iramp,xramp,yramp] = geoimread('ramp2_dem_osu91a200m.tif'); % The OSU version of the file is referenced to the geoid whereas the RAM2_DEM.tif is referenced to wgs84. Use OSU.  
Iramp = double(Iramp); 
Iramp(Iramp==0) = NaN; % Convert zeros to NaN so interpolation won't produce thin ice shelf edges.  
h_ramp = interp2(xramp,yramp,Iramp,X_bm,Y_bm); 
h_ramp(bwdist(~isfinite(h_ramp))*.5<15) = NaN;  % Eliminates 15 km perimeter b/c drooping ice sheet edges in Bamber and RAMP dems

%% Remove firn air content 
% (BedMachine already has FAC removed) 

xf = h5read('FULL_CUBE_v4.h5','/x');
yf = h5read('FULL_CUBE_v4.h5','/y');
fac = interp2(xf,yf,permute(squeeze(mean(h5read('FULL_CUBE_v4.h5','/fac_gemb8'),'omitnan')),[2 1]),X_bm,Y_bm); 

hydro = 9.3364; % factor for converting surface elevation to thickness. 
H_rema = (h_rema - fac)*hydro; 
H_b2 = (h_b2 - fac)*hydro; 
H_bam = (h_bam-fac)*hydro; 
H_ramp = (h_ramp-fac)*hydro; 

% Prevent new grounding: 
H_max = -1.12*bedmachine_data('bed'); % maximum thickness that could exist at neutral buoyancy (no new grounding). if hydrostatic factor thickness2freeboard(1)=0.1071, then H_max = bed/(0.1071-1) = 1.1199*bed.   
H_rema(H_rema>H_max) = H_max(H_rema>H_max); 
H_b2(H_b2>H_max) = H_max(H_b2>H_max); 
H_bam(H_bam>H_max) = H_max(H_bam>H_max); 
H_ramp(H_ramp>H_max) = H_max(H_ramp>H_max); 

% Eliminate negative thickness: 
H_rema(H_rema<0) = nan; 
H_b2(H_b2<0) = nan; 
H_bam(H_bam<0) = nan; 
H_ramp(H_ramp<0) = nan; 

clear h* Iramp fac H_max

%% Combine surface heights 

% Distance in kilometers to BedMachine ocean data: 
D = bwdist(mask_bm==0)*0.5; 
D_thresh = 10; % replace thin ice within 10 km of ocean (for example, where Mertz Glacier tongue broke off and is now thin on its end in the BedMachinne data, overwrite that thin ice with the full thickness from Bedmap2.) 

% Start with BedMachine as the basis: 
source = zeros(size(H_bm),'uint8'); 
source(H_bm>0 | grounded) = 1; 

% Start overwriting BedMachine: 
source(H_bm<(H_rema/2) & D<D_thresh & ~grounded) = 2; 
H_bm(source==2) = H_rema(source==2); 

source(H_bm<(H_b2/2) & D<D_thresh & ~grounded) = 3; 
H_bm(source==3) = H_b2(source==3); 

source(H_bm<(H_bam/2) & D<D_thresh & ~grounded) = 4; 
H_bm(source==4) = H_bam(source==4); 

source(H_bm<(H_ramp/2) & D<D_thresh & ~grounded) = 5; 
H_bm(source==5) = H_ramp(source==5); 

H_bm(H_bm==0 & ~grounded) = nan; 

L = bwlabel(isnan(H_bm)); 
H_bm = regionfill(H_bm,L>1); 
source(L>1) = 6; % interpolated

H = interp2(x_bm,y_bm,H_bm,X,Y); % Interpolate compiled thickness grid
isn = isnan(H); 
H(isn) = interp2(x_bm,y_bm,H_bm,X(isn),Y(isn),'nearest'); % Do a second pass with nearest-neighbor interpolation to get the edgy bits.

H_source = interp2(x_bm,y_bm,source,X,Y,'nearest'); 

clear ax bad D_thresh h H_bm grounded D h* H_bm H_b2 H_rema H_ramp H_bam Iramp source tidal X_bm Y_bm x_bm y_bm xramp yramp fac mask_bm isn xf yf mask_bm

readme = 'created by flow_dem_extend.m. H_source 1=BedMachine v2, 2=REMA, 3=Bedmap2, 4=bamber, 5=ramp2_dem_osu91a200m ';
if save_everything
   save('flow_dem_extend_1.mat','-v7.3') % takes 10 minutes to get here
   disp 'saved 1'
else

   %load('flow_dem_extend_1.mat')
   
   figure
   subplot(1,2,1)
   h = imagescn(x,y,H);
   bedmachine
   ax = gca;
   caxis([0 500])
   subplot(1,2,2)
   h(2) = imagescn(x,y,H_source);
   ax(2) = gca;
   bedmachine
   linkaxes(ax,'xy');
end
 
%% Ice shelf ID mask  
% The mask we're loading was generated by iceshelf_mask_generator.m.

load('iceshelf_mask.mat','iceshelves_dilated')

%% Downscale and Inpaint the flow direction grids. 
% Inpainting H*v (rather than just v) effectively weights the flow directions 
% by how much ice is actually flowing, therefore preventing very thin ice
% from having too much influence on the hypothetical flow of a thick
% glacier. 
% 
% This code section takes about two minutes to inpaint at 1/4 scale. 

sc = 1/4; % scale for resizing velocity (for flow *directions* only) 
[vx_r,x_r,y_r] = demresize(vx,x,y,sc); 
vy_r = imresize(vy,sc); 
H_r = imresize(H,sc); 

Hvx_r = inpaint_nans(H_r.*vx_r,4); 
Hvy_r = inpaint_nans(H_r.*vy_r,4); 
Hv_r = hypot(Hvx_r,Hvy_r); 

% Displacement vectors will be used for the flow directions when extrapolating terminus speeds. 
dx = interp2(x_r,y_r,Hvx_r./Hv_r,X,Y); 
dy = interp2(x_r,y_r,Hvy_r./Hv_r,X,Y); 

figure
imagescn(x_r,y_r,hypot(Hvx_r,Hvy_r))
hold on
bedmachine
q = quiversc(x_r,y_r,Hvx_r,Hvy_r,'k','density',300); 
q.AutoScaleFactor = 3;

%clear H_r vx_r vy_r 

%% Create streamlines for velocity (Round 1) 
% This coarse resolution round is only to constrain the faraway bits later
% when the entire grid is inpainted. There's no real practical scientific
% purpose to this coarse resolution, but it does ultimately help us create
% a beautiful velocity map, and so it's well worth the extra ten minutes or
% so that it takes for this section to run. 
% 
% We're doing multiple rounds simply because of memory and speed
% limitations on my laptop. Ideally we'd do like a million steps in increments 
% of 0.01 or something, but my lil laptop just can't handle it, so this'll
% have to do. 

% Fill any nan holes in the velocity field, but don't fill the ocean: 
L = bwlabel(isnan(vx)); % Label the nan regions, and the ocean will be L=1.  
vx = regionfill(vx,L>1); 
vy = regionfill(vy,L>1);
v_source(L>1) = 4; % interpolated
v = hypot(vx,vy); 

% Use the perimeter of the velocity measurements as the seed locations: 
perim = bwperim(isfinite(vx));
xseed = X(perim); 
yseed = Y(perim); 
[row,col] = find(perim); 
XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.4 6000]);

% Loop through each seed location to extrapolate terminus velocity and ice shelf name along flowlines:  
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v(row(k),col(k)); 
   V{k}(:,2) = iceshelves_dilated(row(k),col(k)); 
end

% Concatenate the cells into nan-separated arrays: 
M = cell2nancat(XY); % cell2nancat is a function at the bottom of this script. 
V = cell2nancat(V); 

% Optional plot to verify things are looking good: 
% skip = 500; 
% figure
% fastscatter(M(1:skip:end,1),M(1:skip:end,2),V(1:skip:end,1),'markersize',1); 
% bedmachine 
% caxis([0 2000]) 

% Are the seed locations captured in M? I'm not sure, so for good measure we'll create these arrays:  
xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
vv = [v(perim);V(:,1)] ;
vvs = [iceshelves_dilated(perim);V(:,2)] ;

% Remove NaNs: 
isf = isfinite(vv);
xx = xx(isf);
yy = yy(isf);
vv = vv(isf);
vvs = vvs(isf);

% Grid up the streamline data:  
vg = gridbin(xx,yy,vv,x,y); % gridbin is on my GitHub. It bins the average velocity observation in each grid cell. It's not as beautiful as gridfit, but it's the only way to deal with nearly a billion scattered data points.   
vsg = gridbin(xx,yy,vvs,x,y,@median); % should really be @mode, but @median is waaay faster 

% Create temporary grids so we don't accidentally do anything dumb and
% overwrite any good data: 
tmpvx = vx; 
tmpvy = vy; 
tmpvs = iceshelves_dilated; 

% Wherever gridbin gives us data, overwrite tmpvx, tmpvy, and tmpvs:  
isf = isfinite(vg);
tmpvx(isf) = dx(isf).*vg(isf); 
tmpvy(isf) = dy(isf).*vg(isf); 
isf = vsg>0; 
tmpvs(isf) = vsg(isf); 

disp 'done first round vel' 

%% Velocity streamlines (Round 2): 
% Slightly denser grid (0.1 instead of 0.4 spacing) does not extend as far
% beyond the current coastlines as Round 1, and doesn't completely fill the ocean
% domain, but does cover up to about 600 km away from the coast, and that's way 
% more than we actually need for our science. 

XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.1 6000]);
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v(row(k),col(k)); 
   V{k}(:,2) = iceshelves_dilated(row(k),col(k)); 
end

M = cell2nancat(XY);
V = cell2nancat(V); 

xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
vv = [v(perim);V(:,1)] ;
vvs = [iceshelves_dilated(perim);V(:,2)] ;
isf = isfinite(vv);
xx = xx(isf);
yy = yy(isf);
vv = vv(isf);
vvs = vvs(isf);

vg = gridbin(xx,yy,vv,x,y);
vsg = gridbin(xx,yy,vvs,x,y,@median);

% Overwrite anythinng that was written by the previous coarse solution:
isf = isfinite(vg);
tmpvx(isf) = dx(isf).*vg(isf); 
tmpvy(isf) = dy(isf).*vg(isf); 
isf = vsg>0;
tmpvs(isf) = vsg(isf); 

disp 'done second round vel'  

%% Velocity streamlines (Round 3): 
% This round is the most dense, and extends about 120 km from the coast. 

XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.02 6000]);
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v(row(k),col(k)); 
   V{k}(:,2) = iceshelves_dilated(row(k),col(k)); 
end

M = cell2nancat(XY);
V = cell2nancat(V); 

xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
vv = [v(perim);V(:,1)] ;
vvs = [iceshelves_dilated(perim);V(:,2)] ;
isf = isfinite(vv);
xx = xx(isf);
yy = yy(isf);
vv = vv(isf);
vvs = vvs(isf);

vg = gridbin(xx,yy,vv,x,y);
vsg = gridbin(xx,yy,vvs,x,y,@median);

isf = isfinite(vg);
tmpvx(isf) = dx(isf).*vg(isf); 
tmpvy(isf) = dy(isf).*vg(isf); 
isf = vsg>0; 
tmpvs(isf) = vsg(isf); 

imagesc(x,y,tmpvs)
axis image off
bedmachine
axis xy

disp 'done third round vel'  

%% Fill remaining bits of missing data: 

% Fill holes in each ice shelf name mask (because there wouldn't be another ice shelf within an ice shelf): 
for k = 1:181
   tmp = imfill(tmpvs==k,'holes'); 
   tmpvs(tmp) = k; 
end

% There are still some random holes in the ice shelf name mask. Fill the little ones with a median filter, and there still will be some holes that sort of look like a zipper, but don't worry about them, it's fine, just deal with it. 
tmpvsmd = medfilt2(tmpvs,[1 1]*5); 
tmpvs(tmpvs==0) = tmpvsmd(tmpvs==0); 

% Fill holes in the velocity data: 
tmpvx2 = regionfill(tmpvx,isnan(tmpvx)); 
tmpvy2 = regionfill(tmpvy,isnan(tmpvy)); 

vx = tmpvx2; 
vy = tmpvy2; 
iceshelves = tmpvs;
iceshelves(bwdist(isfinite(H))*diff(x(1:2)/1000)>100) = 0;  % trims anything more than 100 km from original compiled thickness data.  
v_source(v_source==0 & isfinite(vx)) = 5; 

if false
figure
imagescn(x,y,hypot(vx,vy))
cmocean thermal
hgl = bedmachine('gl','color','k');
hold on
axis off
set(gca,'colorscale','log')
caxis([1.5 4000])
if save_everything
   export_fig flow_extruded3.png -r600
   close all
end
end

%% 

if save_everything
 save('flow_dem_extend_2.mat','vx','vy','iceshelves','v_source','x','y','-v7.3')
end

clear tmp* V vsg vv vvs iceshelves_* bad ch row col D* dx dy L M perim q xx yy

%% Extrude Thickness Round 1 

perim = bwperim(imfill(H>10,'holes'));
xseed = X(perim); 
yseed = Y(perim); 
[row,col] = find(perim); 

Hf = filt2(H,240,5e3,'lp'); 

XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.4 6000]);
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = Hf(row(k),col(k)); 
end

M = cell2nancat(XY);
V = cell2nancat(V); 

xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
hh = [H(perim);V(:,1)] ;
isf = isfinite(hh);
xx = xx(isf);
yy = yy(isf);
hh = hh(isf);

hg = gridbin(xx,yy,hh,x,y);

tmpH = H; 
isf = isfinite(hg);
tmpH(isf) = hg(isf); 

disp 'done thickness round 1' 

%% Extrude Thickness Round 2

XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.1 6000]);
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = Hf(row(k),col(k)); 
end

M = cell2nancat(XY);
V = cell2nancat(V); 

xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
hh = [H(perim);V(:,1)] ;
isf = isfinite(hh);
xx = xx(isf);
yy = yy(isf);
hh = hh(isf);

hg = gridbin(xx,yy,hh,x,y);

isf = isfinite(hg);
tmpH(isf) = hg(isf); 

disp 'done thickness round 2' 

%% Extrude Thickness Round 3

XY = stream2(x_r,y_r,Hvx_r,Hvy_r,xseed,yseed,[0.02 6000]);
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = Hf(row(k),col(k)); 
end

M = cell2nancat(XY);
V = cell2nancat(V); 

xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
hh = [H(perim);V(:,1)] ;
isf = isfinite(hh);
xx = xx(isf);
yy = yy(isf);
hh = hh(isf);

hg = gridbin(xx,yy,hh,x,y);

isf = isfinite(hg);
tmpH(isf) = hg(isf); 

disp 'done thickness round 3' 

%% 

H = regionfill(tmpH,isnan(tmpH)); 
H_source(H_source==0) = 7; 
v_source(v_source==0) = 5; 

%% Write 

if save_everything
   
   % Switch the dimension order to match convention: 
   iceshelves = ipermute(iceshelves,[2 1]); 
   H = ipermute(H,[2 1]); 
   H_source = ipermute(H_source,[2 1]); 
   vx = ipermute(vx,[2 1]); 
   vy = ipermute(vy,[2 1]); 
   v_source = ipermute(v_source,[2 1]); 


   fn = ['extruded_antarctica_',datestr(now,'yyyy-mm-dd'),'.h5'];

   h5create(fn,'/x',size(x),'Datatype','single')
   h5write(fn,'/x',single(x))

   h5create(fn,'/y',size(y),'Datatype','single')
   h5write(fn,'/y',single(y))

   h5create(fn,'/thickness',size(iceshelves),'Datatype','single')
   h5write(fn,'/thickness',single(H))
   h5writeatt(fn,'/thickness','units','meters')
   h5writeatt(fn,'/thickness','Description','Ice thickness compiled from all available datasets and extruded. Hydrostatic assumption (multiplying surface by 9.34) is applied to all floating ice.')

   h5create(fn,'/vx',size(vx),'Datatype','single')
   h5write(fn,'/vx',single(vx))
   h5writeatt(fn,'/vx','units','m/yr')
   h5writeatt(fn,'/vx','Description','Ice velocity in the polar stereographic x direction.') 

   h5create(fn,'/vy',size(vy),'Datatype','single')
   h5write(fn,'/vy',single(vy))
   h5writeatt(fn,'/vy','units','m/yr')
   h5writeatt(fn,'/vy','Description','Ice velocity in the polar stereographic y direction.') 

   h5create(fn,'/thickness_source',size(H_source),'Datatype','uint8')
   h5write(fn,'/thickness_source',H_source)
   h5writeatt(fn,'/thickness_source','data source','1=BedMachine v2, 2=REMA-FAC, 3=Bedmap2-FAC, 4=Bamber-FAC, 5=RAMP2-FAC, 6=interpolated, 7=extrapolated. And FAC is the mean from GEMB.')

   h5create(fn,'/v_source',size(v_source),'Datatype','uint8')
   h5write(fn,'/v_source',v_source)
   h5writeatt(fn,'/v_source','data source','1=ITS_LIVE (Gardner), 2=MEaSUREs v2 (Rignot), 3=error weighted mean of ITS_LIVE and MEaSURES v2, 3=interpolated, 4=extrapolated.')

   h5create(fn,'/iceshelf_mask',size(iceshelves),'Datatype','uint8')
   h5write(fn,'/iceshelf_mask',iceshelves)
   h5writeatt(fn,'/iceshelf_mask','ice shelf names','see iceshelves_2008_v2_names.csv for all 181 ice shelf names.')
   h5writeatt(fn,'/iceshelf_mask','data source','ice shelf masks dilated and extruded from Mouginot 2008 data.')

   h5writeatt(fn,'/','Description','Antarctic ice velocity and thickness extruded along flow at the end of ice shelves. Created by flow_dem_extend.m and on GitHub in a repository named ice-shelf-geometry.')
   h5writeatt(fn,'/','Author','Chad A. Greene, NASA/JPL') 
   h5writeatt(fn,'/','Projection','EPSG:3031 South Polar Stereographic, standard parallel 71S') 

end

%% Extrapolate velocity (old method) 
% 
% L = bwlabel(isnan(vx)); 
% vx = regionfill(vx,L>1); 
% vy = regionfill(vy,L>1);
% 
% vxold = vx; 
% vyold = vy; 
% 
% dx = 500*interp2(x_r,y_r,Hvx_r./Hv_r,X,Y); 
% dy = 500*interp2(x_r,y_r,Hvy_r./Hv_r,X,Y); 
% 
% vx_f = filt2(vx,240,1000,'lp'); 
% vy_f = filt2(vy,240,1000,'lp'); 
% 
% v_source = uint8(isfinite(vx)); 
% 
% iceshelves_dilated(isnan(vx_f)) = 0; 
% st = strel('disk',1);
% for k = 1:50
%    perim = bwperim(imdilate(isfinite(vx), st));
% 
%    vx(perim) = interp2(x,y,vx,X(perim)-dx(perim),Y(perim)-dy(perim)); 
%    vy(perim) = interp2(x,y,vy,X(perim)-dx(perim),Y(perim)-dy(perim)); 
%    
% %    % Where linear interpolation failed (probably due to presence of nans), try nearest-neighbor:
%    isn = perim & isnan(vx); 
%    vx(isn) = interp2(x,y,vx,X(isn)-dx(isn),Y(isn)-dy(isn),'nearest'); 
%    vy(isn) = interp2(x,y,vy,X(isn)-dx(isn),Y(isn)-dy(isn),'nearest'); 
% 
%    iceshelves_dilated(perim) = interp2(x,y,iceshelves_dilated,X(perim)-dx(perim),Y(perim)-dy(perim),'nearest'); 
% 
%    v_source(perim & isfinite(vx)) = 2; 
%    k
% end
% 
% figure
% imagescn(x,y,hypot(vx,vy))
% caxis([0 2000])
% bedmachine
% 
% 
% % save('flow_dem_extend_2.mat','vx','vy','iceshelves_dilated','v_source','x','y','-v7.3') 
% 
% %% Clear up memory 
% 
% clear vx_f vy_f H_r Hv* vx_r vy_r y_r x_r sc tmpx tmpy 
% 
% %% 
% 
% H_f = filt2(H,240,2500,'lp'); 
% 
% st = strel('disk',1);
% for k = 1:416
%    perim = bwperim(imdilate(isfinite(H) | ~ocean, st));
% 
%    tmp = interp2(x,y,H_f,X(perim)-dx(perim),Y(perim)-dy(perim)); 
%    H(perim) = tmp; 
%    H_f(perim) = tmp; 
% 
%    H_source(perim & isfinite(H)) = 2; 
%    k
% end