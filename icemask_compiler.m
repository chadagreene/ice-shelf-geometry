% This script gets ice masks from three different datasets to work 
% nicely together. Run this script *after* running 
% * icemask_modis.m, 
% * icemask_s1a.m, 
% * icemask_ramp2.m 
% 
% Chad Greene, NASA JPL, September 2021. 
%% Load data cubes: 

R = load('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_ramp2b.mat'); 
F = load('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_modis.mat'); 
MOA = load('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_moa.mat'); 
S = load('/Users/cgreene/Documents/MATLAB/DEM_generation/icemask_s1a.mat'); 

D = load('iceshelves_2008_v2.mat'); 

x = F.x; 
y = F.y; 
[X,Y] = meshgrid(x,y); 

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-09-02.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
H = permute(h5read(fn,'/thickness'),[2 1]); 
%H_source = permute(h5read(fn,'/thickness_source'),[2 1]); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 
%v = hypot(vx,vy); 

%%

if true
figure
plot_calving_ts(R,F,MOA,S)
title 'unmodified data'
% export_fig unmodified_iceshelf_area_timeseries.png -r500 -painters
end
%% Remove islands

C = shaperead('scripps_antarctica_polygons_v1.shp'); 

ind = [C.Id]==3; % isolated islands 
islands = inpolygon_map(X,Y,[C(ind).X],[C(ind).Y]); 

% This grounded iceberg by West Ice Shelf shows up in some datasets and not others: 
xi = [2473952,2475928,2478624,2478983,2479882,2480421,2482398,2484914,2487610,2489047,2489766,2490665,2490125,2490485,2490665,2491383,2490844,2491563,2493899,2496595,2499291,2503064,2507198,2508096,2524989,2567759,2600286,2604958,2578721,2540803,2479163,2459036,2464248,2471256]; 
yi = [361600,359084,356748,354591,353154,351356,350278,348122,345786,343090,341473,338777,337519,334644,331948,329971,328354,324940,321345,315954,310204,303375,296546,290256,289537,291155,295288,404011,411918,411379,407784,403112,380828,365913];
berg = inpolygon_map(X,Y,xi,yi); 

grounded = ismember(bedmachine_interp('mask',X,Y),[1 2 4]); 

%%

always_ice = all(cat(3,R.ice,MOA.ice,S.ice,F.ice),3); 

% Mask out any islands or icebergs, but make sure we keep anything tthat's ice in ALL datasets, then put islands back into all the datasets : 
R = mask_out_ice(R,islands |  berg, always_ice, islands | grounded); 
S = mask_out_ice(S,islands | berg, always_ice, islands | grounded); 
MOA = mask_out_ice(MOA,islands  | berg, always_ice, islands | grounded); 
F = mask_out_ice(F,islands | berg, always_ice, islands | grounded); 

%never_ice = ~any(cat(3,R.ice,MOA.ice,S.ice,F.ice),3); 
clear berg 
%%

if true
   figure
   hp = plot_calving_ts(R,F,MOA,S); 
   %export_fig slightly_modified_iceshelf_area_timeseries.png -r500 -painters
end

%%

fraser_only = all(cat(3,F.ice,S.ice),3) & ~any(cat(3,MOA.ice,R.ice),3); 
Fonly = all(F.ice,3) & ~any(cat(3,R.ice,MOA.ice,S.ice),3); 
%Sonly = all(S.ice,3) & ~any(cat(3,R.ice,MOA.ice,F.ice),3); 
%Ronly = all(R.ice,3) & ~any(cat(3,S.ice,MOA.ice,F.ice),3); 
MOAonly = all(MOA.ice,3) & ~any(cat(3,S.ice,R.ice,F.ice),3); 

fraser_miss = ~any(cat(3,F.ice,S.ice),3) & all(cat(3,MOA.ice,R.ice),3); 
Fmiss = ~any(F.ice,3) & all(cat(3,R.ice,MOA.ice,S.ice),3); 
%Smiss = ~any(S.ice,3) & all(cat(3,R.ice,MOA.ice,F.ice),3); 
Rmiss = ~any(R.ice,3) & all(cat(3,S.ice,MOA.ice,F.ice),3); 
MOAmiss = ~any(MOA.ice,3) & all(cat(3,S.ice,R.ice,F.ice),3); 

%% 
% 
% Remove grid cells from F that ONLY exist in F, and add the ones that are missing only from F  
F = mask_out_ice(F,Fonly | fraser_only, Fmiss | fraser_miss, islands | grounded); 
S = mask_out_ice(S,fraser_only,fraser_miss, islands | grounded); 
MOA = mask_out_ice(MOA,MOAonly,MOAmiss, islands | grounded); 
R = mask_out_ice(R,false(size(X)),Rmiss, islands | grounded); 

clear Fonly fraser_only Fmiss fraser_miss MOAonly MOAmiss insignificant never_ice H_source  Rmiss Ronly
%%  
FnotMOA = all(F.ice(:,:,[5 10 15]),3) & ~any(MOA.ice,3); 
MOAnotF = ~any(F.ice(:,:,[5 10 15]),3) & all(MOA.ice,3); 

F = mask_out_ice(F,FnotMOA,MOAnotF, islands | grounded); 

clear FnotMOA MOAnotF
%%

SnotF = all(S.ice(:,:,[1 2 3 7]),3) & ~any(F.ice(:,:,16:19),3); 
FnotS = ~any(S.ice(:,:,[1 2 3 7]),3) & all(F.ice(:,:,16:19),3); 

% consider doing SnotF & v<400

% Most SnotF appears to be digitization/interpretation error in S, so mask it ouut of S. 
% whereas most FnotS is a glob near thwaites and apparently an island near
% the peninsula, so mask it out of F. 

S = mask_out_ice(S,SnotF, always_ice, islands | grounded); 
F = mask_out_ice(F,FnotS, always_ice, islands | grounded); 

%% 
% If we do MOAnotF again, we find that it's mostly edge bits that can be added to every dataset 

MOAnotF = ~any(F.ice(:,:,[5 10 15]),3) & all(MOA.ice,3); 

F = mask_out_ice(F,false(size(mask)),MOAnotF, islands | grounded); 
S = mask_out_ice(S,false(size(mask)),MOAnotF, islands | grounded); 
R = mask_out_ice(R,false(size(mask)),MOAnotF, islands | grounded); 

readme = 'Ice masks adjusted on an entire-campaign basis (not each timestep individually) to make the different mappings roughly consistent with each other. Saved by icemask_compiler.m'; 
%save('/Users/cgreene/Documents/MATLAB/DEM_generation/icemasks_adjusted.mat','R','F','MOA','S','readme','-v7.3')

% figure
% plot_calving_ts(R,F,MOA,S)
clear SnotF MOAnotF C FnotS

%% 

dt = 0.5; % dt is 0.5 year
vscale = 1.1; % Factor of safety for carving and filling. 1 assumes ITS_LIVE velocity measuremets are perfect. 1.1 allows 10% error in the velocities before carving away willy nilly. 
Xb = repmat(single(X),[1 1 10]); % backward-looking 
Yb = repmat(single(Y),[1 1 10]);
Xf = Xb; 
Yf = Yb; 
Xbtmp = X; 
Ybtmp = Y; 
Xftmp = X; 
Yftmp = Y; 
for k = 1:10
   % Location of the grid cell k dt's before
   Xbtmp = Xbtmp-interp2(x,y,vx*dt*vscale,Xb(:,:,k),Yb(:,:,k)); 
   Ybtmp = Ybtmp-interp2(x,y,vy*dt*vscale,Xb(:,:,k),Yb(:,:,k)); 
   Xb(:,:,k) = Xbtmp; 
   Yb(:,:,k) = Ybtmp; 
   
   Xftmp = Xftmp+interp2(x,y,vx*dt*vscale,Xf(:,:,k),Yf(:,:,k)); 
   Yftmp = Yftmp+interp2(x,y,vy*dt*vscale,Xf(:,:,k),Yf(:,:,k)); 
   Xf(:,:,k) = Xftmp; 
   Yf(:,:,k) = Yftmp; 
end
clear Ybtmp Xbtmp Xftmp Yftmp vx vy
disp yep 

%% Create ice cube 

year = [1997.75 2000.2 2000.75 2001.2:2021.2]; 

ice = false(size(X,1),size(X,2),length(year)); 

% Start by filling with baseline Radarsat and MODIS mappings: 
ice(:,:,1) = R.ice(:,:,1); 
ice(:,:,[7 12 17]) = MOA.ice; 

figure
plot_calving_ts(R,F,MOA,S,ice,year)

%% Adjust RAMP 2000
% It is clear from inspection that the 1997 coastline is better than the 2000 coastline. It was 
% clearly more carefully made, as the 1997 coastline features detailed
% rifts whereas the 2000 coastline appears to smooth over them. Therefore,
% if any ice suddenly "appears" in the 2000 Ramp coastline that couldn't
% have gotten there after just 3 years of advection, remove it from the
% 2000 mapping. Then we'll be able to confidently use the 2000 mapping as a
% reference for Modis mappings. 

was_ice = interp2(x,y,double(ice(:,:,1)),Xb(:,:,6),Yb(:,:,6))>0.5; 
will_be_ice = interp2(x,y,double(ice(:,:,7)),Xf(:,:,7),Yf(:,:,7))>0.5; 

tmp = R.ice(:,:,2); 
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,3) = tmp; 

%% Adjust MOA

% Adjust MOA2004 based on ramp2000 and MOA2009
will_be_ice = interp2(x,y,double(ice(:,:,12)),Xf(:,:,10),Yf(:,:,10))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,3)),Xb(:,:,7),Yb(:,:,7))>0.5; 

tmp = ice(:,:,7);
tmp(will_be_ice & was_ice) = true;   
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,7) = tmp; 

% Adjust MOA2009 based on adjusted MOA2004 and MOA2014
will_be_ice = interp2(x,y,double(ice(:,:,17)),Xf(:,:,10),Yf(:,:,10))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,7)),Xb(:,:,10),Yb(:,:,10))>0.5; 

tmp = ice(:,:,12);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,12) = tmp; 

% Adjust MOA2014 based on adjusted MOA2009. It's a high bar for adding ice to MOA2014: It must be in the Modis record 3 years later AND the Sentinel record 5 years later AND it must've been ice in the previous MOA.  
was_ice = interp2(x,y,double(ice(:,:,12)),Xb(:,:,10),Yb(:,:,10))>0.5; 
will_be_ice = interp2(x,y,double(S.ice(:,:,5)),Xf(:,:,10),Yf(:,:,10))>0.5 & interp2(x,y,double(F.ice(:,:,18)),Xf(:,:,6),Yf(:,:,6))>0.5; 

tmp = ice(:,:,17);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,17) = tmp; 

%% Do all the Sentinel: 

for k = 1:7
   was_ice = interp2(x,y,double(ice(:,:,16+k)),Xb(:,:,2),Yb(:,:,2))>0.5; 

   switch k
      case 1 % 2015
         % Looking ahead to only add ice where modis will have ice in 2 years AND sentinel will have ice in 5 years:  
         will_be_ice = interp2(x,y,double(S.ice(:,:,6)),Xf(:,:,10),Yf(:,:,10))>0.5 & interp2(x,y,double(F.ice(:,:,18)),Xf(:,:,4),Yf(:,:,4))>0.5; 
      case 2 % 2016
         will_be_ice = interp2(x,y,double(S.ice(:,:,7)),Xf(:,:,10),Yf(:,:,10))>0.5 & interp2(x,y,double(F.ice(:,:,19)),Xf(:,:,10),Yf(:,:,10))>0.5; 
      case 3 % 2017
         will_be_ice = interp2(x,y,double(S.ice(:,:,7)),Xf(:,:,8),Yf(:,:,8))>0.5 & interp2(x,y,double(F.ice(:,:,19)),Xf(:,:,8),Yf(:,:,8))>0.5; 
      case 4 % 2018
         will_be_ice = interp2(x,y,double(S.ice(:,:,7)),Xf(:,:,6),Yf(:,:,6))>0.5 & interp2(x,y,double(F.ice(:,:,19)),Xf(:,:,6),Yf(:,:,6))>0.5; 
      case 5 % 2019
         will_be_ice = interp2(x,y,double(S.ice(:,:,7)),Xf(:,:,4),Yf(:,:,4))>0.5 & interp2(x,y,double(F.ice(:,:,19)),Xf(:,:,4),Yf(:,:,4))>0.5; 
      case 6 % 2020
         will_be_ice = interp2(x,y,double(S.ice(:,:,7)),Xf(:,:,2),Yf(:,:,2))>0.5 & interp2(x,y,double(F.ice(:,:,19)),Xf(:,:,2),Yf(:,:,2))>0.5; 
      case 7
         will_be_ice = false(size(X)); 
   end
   
   tmp = S.ice(:,:,k);
   tmp(will_be_ice & was_ice) = true; 
   tmp(~was_ice & ~will_be_ice) = false; 
   tmp(grounded | islands) = true; 
   tmp = imfill(tmp,8,'holes'); 
   tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
   tmp(grounded | islands) = true; 
   ice(:,:,17+k) = tmp; 
   k
end


s1 = squeeze(sum(sum(ice)))*.24^2;
isf = s1>0;
plot(year(isf),s1(isf),'v-')

%% Adjust the first MODIS based mapping 

% First modis mapping is 2.5 years after RAMP1 and 0.5 years before RAMP2
will_be_ice = interp2(x,y,double(ice(:,:,3)),Xf(:,:,1),Yf(:,:,1))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,1)),Xb(:,:,5),Yb(:,:,5))>0.5; 

tmp = F.ice(:,:,1); 
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,2) = tmp; 

%% Go back and make sure Radarsat 1997 contains the ice seen in modis 2000 and ramp 2000 
% Make sure any ice we add to 1997 is present in *both* 2000 datasets. 

will_be_ice = interp2(x,y,double(ice(:,:,2)),Xf(:,:,5),Yf(:,:,5))>0.5; 
will_be_ice2 = interp2(x,y,double(ice(:,:,3)),Xf(:,:,6),Yf(:,:,6))>0.5; 
tmp = ice(:,:,1); 
tmp(will_be_ice & will_be_ice2) = true; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,1) = tmp; 

%% The rest of the MODIS

% The 2001 modis mapping was 6 months after ramp2 and 3 years before 2004
% moa mapping
will_be_ice = interp2(x,y,double(ice(:,:,7)),Xf(:,:,6),Yf(:,:,6))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,3)),Xb(:,:,1),Yb(:,:,1))>0.5; 
tmp = F.ice(:,:,2);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,4) = tmp; 

% Modis 2002
will_be_ice = interp2(x,y,double(ice(:,:,7)),Xf(:,:,4),Yf(:,:,4))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,4)),Xb(:,:,2),Yb(:,:,2))>0.5; 
tmp = F.ice(:,:,3);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,5) = tmp; 

% Modis 2003
will_be_ice = interp2(x,y,double(ice(:,:,7)),Xf(:,:,2),Yf(:,:,2))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,5)),Xb(:,:,2),Yb(:,:,2))>0.5; 
tmp = F.ice(:,:,4);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,6) = tmp; 


% 2005: 
will_be_ice = interp2(x,y,double(ice(:,:,12)),Xf(:,:,8),Yf(:,:,8))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,7)),Xb(:,:,2),Yb(:,:,2))>0.5; 
tmp = F.ice(:,:,6);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,8) = tmp; 

% 2006:
will_be_ice = interp2(x,y,double(ice(:,:,12)),Xf(:,:,6),Yf(:,:,6))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,7)),Xb(:,:,4),Yb(:,:,4))>0.5; 
tmp = F.ice(:,:,7);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,9) = tmp; 

% 2007:
will_be_ice = interp2(x,y,double(ice(:,:,12)),Xf(:,:,4),Yf(:,:,4))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,7)),Xb(:,:,6),Yb(:,:,6))>0.5; 
tmp = F.ice(:,:,8);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,10) = tmp; 

% 2008:
will_be_ice = interp2(x,y,double(ice(:,:,12)),Xf(:,:,2),Yf(:,:,2))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,7)),Xb(:,:,8),Yb(:,:,8))>0.5; 
tmp = F.ice(:,:,9);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,11) = tmp; 

% 2010:
will_be_ice = interp2(x,y,double(ice(:,:,17)),Xf(:,:,8),Yf(:,:,8))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,12)),Xb(:,:,2),Yb(:,:,2))>0.5; 
tmp = F.ice(:,:,11);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,13) = tmp; 

% 2011:
will_be_ice = interp2(x,y,double(ice(:,:,17)),Xf(:,:,6),Yf(:,:,6))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,12)),Xb(:,:,4),Yb(:,:,4))>0.5; 
tmp = F.ice(:,:,12);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,14) = tmp; 

% 2012: 
will_be_ice = interp2(x,y,double(ice(:,:,17)),Xf(:,:,4),Yf(:,:,4))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,12)),Xb(:,:,6),Yb(:,:,6))>0.5; 
tmp = F.ice(:,:,13);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,15) = tmp; 

% 2013: 
will_be_ice = interp2(x,y,double(ice(:,:,17)),Xf(:,:,2),Yf(:,:,2))>0.5; 
was_ice = interp2(x,y,double(ice(:,:,12)),Xb(:,:,8),Yb(:,:,8))>0.5; 
tmp = F.ice(:,:,14);
tmp(will_be_ice & was_ice) = true; 
tmp(~was_ice & ~will_be_ice) = false; 
tmp(grounded | islands) = true; 
tmp = imfill(tmp,8,'holes'); 
tmp(~bwselect(tmp,floor(length(x)/2),floor(length(y)/2))) = false; 
tmp(grounded | islands) = true; 
ice(:,:,16) = tmp; 


%%
figure
plot_calving_ts(R,F,MOA,S,ice,year)

%% Re-ground any islands: 

% grounded = isgrounded(X,Y);
% 
% for k = 1:length(year)
%    tmp = ice(:,:,k); 
%    tmp(grounded | islands) = true; 
%    ice(:,:,k) = tmp; 
% end


for k = 1:length(year)
   
   [cx{k},cy{k}] = mask2outline(x,y,ice(:,:,k)); 
k
end

readme = 'Adjusted ice grid created by ice_compiler.m'; 
%save('icemask_composite.mat','ice','x','y','year','readme','cx','cy','-v7.3')

return

%% Write 

if save_everything
   
   % Switch the dimension order to match convention: 
   ice = ipermute(ice,[2 1 3]); 

   fn = ['/Users/cgreene/Documents/MATLAB/DEM_generation/itslive_coastmaskb_',datestr(now,'yyyy-mm-dd'),'.h5'];

   h5create(fn,'/x',size(x),'Datatype','single')
   h5write(fn,'/x',single(x))

   h5create(fn,'/y',size(y),'Datatype','single')
   h5write(fn,'/y',single(y))

   h5create(fn,'/year',size(year),'Datatype','single')
   h5write(fn,'/year',single(year))
   
   h5create(fn,'/ice',size(ice),'Datatype','int8')
   h5write(fn,'/ice',int8(ice))
   h5writeatt(fn,'/ice','units','boolean')
   h5writeatt(fn,'/ice','Description','True where ice or continent was observed. False otherwise. Island extents remain constant and do not have evolving coastlines.')

   h5writeatt(fn,'/','Description','A physically consistent composite of coastlines from Radarsat, MOA, MODIS (Fraser), and Sentinel 1a. Created by icemask_compiler.m and on GitHub in a repository named ice-shelf-geometry.')
   h5writeatt(fn,'/','Author','Chad A. Greene, NASA/JPL') 
   h5writeatt(fn,'/','Projection','EPSG:3031 South Polar Stereographic, standard parallel 71S') 

end

%%
k=1;

figure
h = imagescn(x,y,ice(:,:,k)); 
axis tight equal off 
col = parula(3); 
cmocean ice 
bedmachine('gl','color',.5*[1 1 1])
modismoaps('coast','color',col(1,:),'year',2004)
modismoaps('coast','color',col(2,:),'year',2009)
modismoaps('coast','color',col(3,:),'year',2014)
title(num2str(floor(year(k))))

%% 
k=k+1; 
for k=1:24
h.CData = ice(:,:,k); 
title(num2str(floor(year(k))))
drawnow
end

%%

always_ice = all(cat(3,R.ice,F.ice,S.ice,MOA.ice,ice),3); 
never_ice = ~any(cat(3,R.ice,F.ice,S.ice,MOA.ice,ice),3); 
%%

iceshelf = 1; 

msk = mask==iceshelf & ~always_ice & ~never_ice; 

AR = local(R.ice,msk,@sum)*.24^2; 
AMOA = local(MOA.ice,msk,@sum)*.24^2; 
AS = local(S.ice,msk,@sum)*.24^2; 
AF = local(F.ice,msk,@sum)*.24^2; 
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

export_fig(['iceshelf_area_timeseries_',D.name{iceshelf},'.png'],'-r300','-painters')
return

%%  Subfunctions 

function S = mask_out_ice(S,mask,keep,islands) 

assert(nargin==4,'Not enough inputs.')

   [m,n,o] = size(S.ice); 
   for k = 1:o
      tmp = S.ice(:,:,k); 
      tmp(mask) = false; 
      if nargin>2
         tmp(keep) = true; 
      end
      tmp = imfill(tmp,8,'holes'); 
      tmp(~bwselect(tmp,floor(n/2),floor(m/2))) = false; 
      if nargin>3
         tmp(islands) = true; 
      end
      S.ice(:,:,k) = tmp; 
   end
end


function h = plot_calving_ts(R,F,MOA,S,ice,year)

   r1 = squeeze(sum(sum(R.ice)))*0.24^2;
   f1 = squeeze(sum(sum(F.ice)))*0.24^2;
   moa1 = squeeze(sum(sum(MOA.ice)))*0.24^2;
   s1 = squeeze(sum(sum(S.ice)))*0.24^2;
   if nargin>4
      ice1 = squeeze(sum(sum(ice)))*0.24^2;
   end
   

   h(1) = plot(R.yr+.75,r1,'o-','linewidth',1);
   hold on
   h(2) = plot(MOA.yr+.2,moa1,'o-','linewidth',1);
   h(3) = plot(S.yr+.2,s1,'o-','linewidth',1);
   h(4) = plot(F.years+.2,f1,'o-','linewidth',1);
   
   switch nargin 
      case 4
         legend('Radarsat','MOA','Fraser S1a','Fraser MODIS','location','best')
      case 6 
         isf = ice1>0; 
         h(5) = plot(year(isf),ice1(isf),'kx-','linewidth',2);
         legend('Radarsat','MOA','Fraser S1a','Fraser MODIS','composite','location','best')
      otherwise
         error('wrong number of inputs') 
   end
   
   box off 
   axis tight
   set(gcf,'renderer','painters') 
   ylabel 'total ice area (km^2)' 
   
   if nargout==0
      clear h
   end

end