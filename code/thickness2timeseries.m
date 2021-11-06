
%% SECOND TRY 


% This script creates thickness_and_melt_rates.mat. 
% I know the filename of this script says timeseries, but for now it's just
% mean rates for 1997-2017. 
% 
% This script takes 4 or 5 min to run. 
% 
% Chad Greene, November 2021. 

%% Load data 

load('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','ground')
load('icemask_composite.mat')
always_ice = imfill(all(ice(:,:,year<=2018),3),8,'holes'); % Take the mean ice mask over the observed period. Taking the sum is faster than mean.  
clear ice year

fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 

%[X,Y] = meshgrid(x,y); 
rho_ice = 917;
D = load('iceshelves_2008_v2.mat');
%tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]) & ~ground; 

%%

yrf = h5read('FULL_CUBE_v4.h5','/t'); 
xf = h5read('FULL_CUBE_v4.h5','/x');
yf = h5read('FULL_CUBE_v4.h5','/y');
[Xf,Yf] = meshgrid(xf,yf); 

Hf = permute(h5read('FULL_CUBE_v4.h5','/H_filt10'),[3 2 1]); 
Hf_err = permute(h5read('FULL_CUBE_v4.h5','/H_err10'),[3 2 1]); 
melt = permute(h5read('FULL_CUBE_v4.h5','/dHdt_melt10'),[3 2 1]); 
melt_errf = permute(h5read('FULL_CUBE_v4.h5','/dHdt_melt_err10'),[3 2 1]); 

smb = permute(h5read('FULL_CUBE_v4.h5','/smb_gemb8'),[3 2 1]); 
smb_errf = permute(h5read('FULL_CUBE_v4.h5','/smb_gemb_err8'),[3 2 1]); 

mask = interp2(x,y,mask,Xf,Yf,'nearest'); 
always_ice = interp2(x,y,double(always_ice),Xf,Yf,'nearest')==1; 

%% 

melt_mean = mean(melt(:,:,yrf>=1997.75),3,'omitnan'); 
melt_err = mean(melt_errf(:,:,yrf>=1997.75),3,'omitnan'); 

smb_mean = mean(smb(:,:,yrf>=1997.75),3,'omitnan'); 
smb_err = mean(smb_errf(:,:,yrf>=1997.75),3,'omitnan'); 

H_trend = trend(Hf(:,:,yrf>=1997.75),yrf(yrf>=1997.75)); 

% Equation 8.17 on pg 188 of Introduction to Error Analysis (Taylor)
del = sum(yrf>=1997.75)*sum(yrf(yrf>=1997.75).^2) - sum(yrf(yrf>=1997.75)).^2; 

H_trend_err = mean(Hf_err(:,:,yrf>=1997.75),3,'omitnan').*sqrt(sum(yrf>=1997.75)/del);

clear melt melt_errf smb smb_errf Hf Hf_err

%% 

valid_data = isfinite(melt_mean) & always_ice; 

% Determine the mass of ice in each grid cell: 
[Lat,~] = ps2ll(Xf,Yf); 
F = rho_ice*1e-12*(median(diff(xf(1:2))) ./ psdistortion(Lat) ).^2; % F is a factor that converts vertical rates of ice change to Gt per grid cell.  
clear Lat

melt_mean2 = cube2rect(melt_mean,valid_data); 
melt_err2 = cube2rect(melt_err,valid_data); 
smb_mean2 = cube2rect(smb_mean,valid_data); 
smb_err2 = cube2rect(smb_err,valid_data); 
H_trend2 = cube2rect(H_trend,valid_data); 
H_trend_err2 = cube2rect(H_trend_err,valid_data); 
F2 = cube2rect(F,valid_data); 
mask2 = cube2rect(mask,valid_data); 
has_data2 = cube2rect(valid_data,valid_data); 

%%

MR = nan(183,1); % meltrate 
MR_err = MR; 
SMB = MR; 
SMB_err = MR; 
dHdt = MR;
dHdt_err = MR; 

for k = 1:181
   if any(mask2==k & has_data2)
      MR(k) = sum(melt_mean2(mask2==k).*F2(mask2==k)); 
      %MR_err(k) = sum(??(mask2==k).*F2(mask2==k)); 
      SMB(k) = sum(smb_mean2(mask2==k).*F2(mask2==k)); 
      SMB_err(k) = sum(smb_err2(mask2==k).*F2(mask2==k)); 
      dHdt(k) = sum(H_trend2(mask2==k).*F2(mask2==k)); 
      dHdt_err(k) = sum(H_trend_err2(mask2==k).*F2(mask2==k)); 
   end
end

MR(182) = sum(melt_mean2(mask2==182).*F2(mask2==182)); 
%MR_err(182) = sum(??(mask2==182).*F2(mask2==182)); 
SMB(182) = sum(smb_mean2(mask2==182).*F2(mask2==182)); 
SMB_err(182) = sum(smb_err2(mask2==182).*F2(mask2==182)); 
dHdt(182) = sum(H_trend2(mask2==182).*F2(mask2==182)); 
dHdt_err(182) = sum(H_trend_err2(mask2==182).*F2(mask2==182)); 

MR(183) = sum(melt_mean2.*F2); 
%MR_err(183) = sum(??.*F2); 
SMB(183) = sum(smb_mean2.*F2); 
SMB_err(183) = sum(smb_err2.*F2); 
dHdt(183) = sum(H_trend2.*F2); 
dHdt_err(183) = sum(H_trend_err2.*F2); 


readme = 'total ice shelf annual rates (Gt/yr) due to MR=meltrate, dHdt=thickness trend, SMB=surface mass balance. Created by thickness2timeseries.m.' 

% save('thickness_and_melt_rates.mat','MR','MR_err','SMB','SMB_err','dHdt','dHdt_err','readme')




%% OLD VERSION 
% INITIALLY I TRIED INTERPOLATING TO THE 240 M GRID, BUT THEN I BEGAN TO
% QUESTION WHETHER THAT'S HELPFUL OR EVEN REASONABLE. ANYWAYS, HERE'S THE
% OLD VERSION: 
% 
% 
% % This script creates thickness_and_melt_rates.mat. 
% % I know the filename of this script says timeseries, but for now it's just
% % mean rates for 1997-2017. 
% % 
% % This script takes 4 or 5 min to run. 
% % 
% % Chad Greene, November 2021. 
% 
% %% Load data 
% 
% load('/Users/cgreene/Documents/MATLAB/DEM_generation/issm_calving_melt_setup.mat','ground')
% load('icemask_composite.mat')
% always_ice = imfill(all(ice(:,:,year<=2018),3),8,'holes'); % Take the mean ice mask over the observed period. Taking the sum is faster than mean.  
% clear ice year
% 
% fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
% mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
% 
% [X,Y] = meshgrid(x,y); 
% rho_ice = 917;
% D = load('iceshelves_2008_v2.mat');
% tidal = ismember(bedmachine_interp('mask',X,Y),[0 3]) & ~ground; 
% 
% %%
% 
% yrf = h5read('FULL_CUBE_v4.h5','/t'); 
% xf = h5read('FULL_CUBE_v4.h5','/x');
% yf = h5read('FULL_CUBE_v4.h5','/y');
% [Xf,Yf] = meshgrid(xf,yf); 
% 
% Hf = permute(h5read('FULL_CUBE_v4.h5','/H_filt10'),[3 2 1]); 
% Hf_err = permute(h5read('FULL_CUBE_v4.h5','/H_err10'),[3 2 1]); 
% melt = permute(h5read('FULL_CUBE_v4.h5','/dHdt_melt10'),[3 2 1]); 
% melt_errf = permute(h5read('FULL_CUBE_v4.h5','/dHdt_melt_err10'),[3 2 1]); 
% 
% smb = permute(h5read('FULL_CUBE_v4.h5','/smb_gemb8'),[3 2 1]); 
% smb_errf = permute(h5read('FULL_CUBE_v4.h5','/smb_gemb_err8'),[3 2 1]); 
% 
% %% 
% 
% melt_mean = mean(melt(:,:,yrf>=1997.75),3,'omitnan'); 
% melt_err = mean(melt_errf(:,:,yrf>=1997.75),3,'omitnan'); 
% 
% smb_mean = mean(smb(:,:,yrf>=1997.75),3,'omitnan'); 
% smb_err = mean(smb_errf(:,:,yrf>=1997.75),3,'omitnan'); 
% 
% H_trend = trend(Hf(:,:,yrf>=1997.75),yrf(yrf>=1997.75)); 
% 
% % Equation 8.17 on pg 188 of Introduction to Error Analysis (Taylor)
% del = sum(yrf>=1997.75)*sum(yrf(yrf>=1997.75).^2) - sum(yrf(yrf>=1997.75)).^2; 
% 
% H_trend_err = mean(Hf_err(:,:,yrf>=1997.75),3,'omitnan').*sqrt(sum(yrf>=1997.75)/del);
% 
% clear melt melt_errf smb smb_errf Hf Hf_err
% 
% %% Fill edges of Fernando's data 
%   
% % Constrain regionfill by setting dH anomalies more than 10 km from data to 0.  
% dist2data_km = 3*bwdist(all(isfinite(H_trend),3));
% 
% melt_mean(dist2data_km>10) = 0; 
% melt_err(dist2data_km>10) = 0; 
% smb_mean(dist2data_km>10) = 0; 
% smb_err(dist2data_km>10) = 0; 
% H_trend(dist2data_km>10) = 0; 
% H_trend_err(dist2data_km>10) = 0; 
% 
% % Fill in the gaps: 
% melt_mean = regionfill(melt_mean,isnan(melt_mean)); 
% melt_err = regionfill(melt_err,isnan(melt_err)); 
% smb_mean = regionfill(smb_mean,isnan(smb_mean)); 
% smb_err = regionfill(smb_err,isnan(smb_err)); 
% H_trend = regionfill(H_trend,isnan(H_trend)); 
% H_trend_err = regionfill(H_trend_err,isnan(H_trend_err)); 
% 
% %% Interpolate to ITS_LIVE 240 m grid
% 
% melt_mean = interp2(xf,yf,melt_mean,X,Y); 
% melt_err = interp2(xf,yf,melt_err,X,Y); 
% smb_mean = interp2(xf,yf,smb_mean,X,Y); 
% smb_err = interp2(xf,yf,smb_err,X,Y); 
% H_trend = interp2(xf,yf,H_trend,X,Y); 
% H_trend_err = interp2(xf,yf,H_trend_err,X,Y); 
% 
% % Even though we're extrapolating a little bit, we'll ultimately only want  
% has_data = interp2(xf,yf,double(dist2data_km==0),X,Y,'nearest'); 
% 
% %%
% 
% % Determine the mass of ice in each grid cell: 
% [Lat,~] = ps2ll(X,Y); 
% F = rho_ice*1e-12*(diff(x(1:2)) ./ psdistortion(Lat) ).^2; % F is a factor that converts vertical rates of ice change to Gt per grid cell.  
% clear Lat
% 
% % A mask of finite ice shelf grid cells:
% isf = tidal & always_ice;
% 
% melt_mean2 = cube2rect(melt_mean,isf); 
% melt_err2 = cube2rect(melt_err,isf); 
% smb_mean2 = cube2rect(smb_mean,isf); 
% smb_err2 = cube2rect(smb_err,isf); 
% H_trend2 = cube2rect(H_trend,isf); 
% H_trend_err2 = cube2rect(H_trend_err,isf); 
% F2 = cube2rect(F,isf); 
% mask2 = cube2rect(mask,isf); 
% has_data2 = cube2rect(has_data,isf); 
% 
% %%
% 
% MR = nan(183,1); % meltrate 
% MR_err = MR; 
% SMB = MR; 
% SMB_err = MR; 
% dHdt = MR;
% dHdt_err = MR; 
% 
% for k = 1:181
%    if any(mask2==k & has_data2)
%       MR(k) = sum(melt_mean2(mask2==k).*F2(mask2==k)); 
%       %MR_err(k) = sum(??(mask2==k).*F2(mask2==k)); 
%       SMB(k) = sum(smb_mean2(mask2==k).*F2(mask2==k)); 
%       SMB_err(k) = sum(smb_err2(mask2==k).*F2(mask2==k)); 
%       dHdt(k) = sum(H_trend2(mask2==k).*F2(mask2==k)); 
%       dHdt_err(k) = sum(H_trend_err2(mask2==k).*F2(mask2==k)); 
%    end
% end
% 
% MR(182) = sum(melt_mean2(mask2==182).*F2(mask2==182)); 
% %MR_err(182) = sum(??(mask2==182).*F2(mask2==182)); 
% SMB(182) = sum(smb_mean2(mask2==182).*F2(mask2==182)); 
% SMB_err(182) = sum(smb_err2(mask2==182).*F2(mask2==182)); 
% dHdt(182) = sum(H_trend2(mask2==182).*F2(mask2==182)); 
% dHdt_err(182) = sum(H_trend_err2(mask2==182).*F2(mask2==182)); 
% 
% MR(183) = sum(melt_mean2.*F2); 
% %MR_err(183) = sum(??.*F2); 
% SMB(183) = sum(smb_mean2.*F2); 
% SMB_err(183) = sum(smb_err2.*F2); 
% dHdt(183) = sum(H_trend2.*F2); 
% dHdt_err(183) = sum(H_trend_err2.*F2); 
% 
% 
% readme = 'total ice shelf annual rates (Gt/yr) due to MR=meltrate, dHdt=thickness trend, SMB=surface mass balance. Created by thickness2timeseries.m.' 
% 
% % save('thickness_and_melt_rates.mat','MR','MR_err','SMB','SMB_err','dHdt','dHdt_err','readme')
% 



