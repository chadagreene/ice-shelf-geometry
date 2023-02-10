% This script just loads and resaves extruded velocity data in a compressed
% NetCDF format. 
% Chad Greene, July 2022. 


fn = '/Users/cgreene/Documents/MATLAB/DEM_generation/extruded_antarctica_2021-10-18.h5'; 
h5disp(fn)

%%

x = h5read(fn,'/x'); 
y = h5read(fn,'/y'); 
vx = permute(h5read(fn,'/vx'),[2 1]); 
vy = permute(h5read(fn,'/vy'),[2 1]); 
v_source = permute(h5read(fn,'/v_source'),[2 1]); 

iceshelf_mask = permute(h5read(fn,'/iceshelf_mask'),[2 1]); 
thickness = permute(h5read(fn,'/thickness'),[2 1]); 
thickness_source = permute(h5read(fn,'/thickness_source'),[2 1]); 

load('icemask_composite.mat'); 
%% Save data 

filename = '/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/antarctica_icemasks_extruded_velocity_and_thickness.nc'; 

disp 'saving data now...'
thickness = uint16(ceil(thickness)); % NC_USHORT
vx = int16(round(vx)); % NC_SHORT
vy = int16(round(vy)); % NC_SHORT
v_source = uint8(v_source); 
thickness_source = uint8(thickness_source); 
iceshelf_mask = uint8(iceshelf_mask); 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
%mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(filename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Antarctic ice masks, extruded velocity and extruded thickness');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','This dataset imagines an Antarctic Ice Sheet without calving fronts. Velocity and thickness observations are extrapolated beyond present-day extents to fill the entire domain. Thickness and velocity extrapolations are constant, and therefore mass rates are somewhat conserved, except that flow may converge or diverge, in which case volume is not conserved. Close to the present-day calving fronts, however, this dataset provides a reasonable first-order approximation of thickness and velocity beyond their observed extents.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date_created',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Institution','NASA Jet Propulsion Laboratory (JPL), California Institute of Technology');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Citation','Chad A. Greene, Alex S. Gardner, Nicole-Jeanne Schlegel & Alexander D. Fraser, 2022. Antarctic calving loss rivals ice-shelf thinning. Nature. https://doi.org/10.1038/s41586-022-05037-w');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GitHub','https://github.com/chadagreene/ice-shelf-geometry');

% 2. Define dimensions
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'inverse_flattening',298.257224);
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
netcdf.putAtt(ncid,mapping_var_id,'scale_factor_at_projection_origin',1);
netcdf.putAtt(ncid,mapping_var_id,'semi_major_axis',6378137);
netcdf.putAtt(ncid,mapping_var_id,'spatial_epsg',3031);
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj','+proj=stere +lat_0=90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs');
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',-71);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',0);

% Define x 
x_id     = netcdf.defDim(ncid,'x',length(x));
x_var_id = netcdf.defVar(ncid,'x','NC_INT',x_id);
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate');
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'units',        'meter');

% Define y
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_INT',y_id);
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate');
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'units',        'meter');

% Define year
year_id     = netcdf.defDim(ncid,'year',length(year));
year_var_id = netcdf.defVar(ncid,'year','NC_FLOAT',year_id);
netcdf.putAtt(ncid,y_var_id,'long_name',    'time');
netcdf.putAtt(ncid,y_var_id,'units',        'decimal year');

% Define vx
vx_var_id = netcdf.defVar(ncid,'vx','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid, vx_var_id,'long_name','projected x component of velocity');
netcdf.putAtt(ncid, vx_var_id,'standard_name','vx');
netcdf.putAtt(ncid, vx_var_id,'units',    'm/yr');
netcdf.putAtt(ncid, vx_var_id,'grid_mapping', 'mapping');

% Define vy
vy_var_id = netcdf.defVar(ncid,'vy','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid, vy_var_id,'long_name','projected y component of velocity');
netcdf.putAtt(ncid, vy_var_id,'standard_name','vy');
netcdf.putAtt(ncid, vy_var_id,'units',    'm/yr');
netcdf.putAtt(ncid, vy_var_id,'grid_mapping', 'mapping');

% Define v_source
v_source_var_id = netcdf.defVar(ncid,'v_source','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid, v_source_var_id,'long_name',    '1=ITS_LIVE v1 (Gardner), 2=MEaSUREs v2 (Rignot), 3=error weighted mean of ITS_LIVE and MEaSURES v2, 3=interpolated, 4=extrapolated.');
netcdf.putAtt(ncid, v_source_var_id,'grid_mapping', 'mapping');

% Define thickness
thick_var_id = netcdf.defVar(ncid,'thickness','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,thick_var_id,'long_name','ice thickness');
netcdf.putAtt(ncid,thick_var_id,'standard_name','land_ice_thickness');
netcdf.putAtt(ncid,thick_var_id,'units',    'meters');
netcdf.putAtt(ncid,thick_var_id,'grid_mapping', 'mapping');

% Define thickness_source
thickness_source_var_id = netcdf.defVar(ncid,'thickness_source','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid, thickness_source_var_id,'long_name',    '1=BedMachine v2, 2=REMA-FAC, 3=Bedmap2-FAC, 4=Bamber-FAC, 5=RAMP2-FAC, 6=interpolated, 7=extrapolated. And FAC is the mean from GEMB.');
netcdf.putAtt(ncid, thickness_source_var_id,'grid_mapping', 'mapping');

% Define ice
ice_var_id = netcdf.defVar(ncid,'ice','NC_BYTE',[x_id y_id year_id]);
netcdf.putAtt(ncid, ice_var_id,'long_name','Binary grid indicates presence of ice.');
netcdf.putAtt(ncid, ice_var_id,'grid_mapping', 'mapping');

% Define mask
mask_var_id = netcdf.defVar(ncid,'iceshelf_mask','NC_UBYTE',[x_id y_id]);
netcdf.putAtt(ncid, mask_var_id,'long_name',    'See greene_Supplementary_Table_1.xlsx for all 181 ice shelf names.');
netcdf.putAtt(ncid, mask_var_id,'data_source',    'Ice shelf masks dilated and extruded from Mouginot 2008 data.');
netcdf.putAtt(ncid, mask_var_id,'grid_mapping', 'mapping');

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,vx_var_id,true,true,9);
netcdf.defVarDeflate(ncid,vy_var_id,true,true,9);
netcdf.defVarDeflate(ncid,v_source_var_id,true,true,9);
netcdf.defVarDeflate(ncid,thick_var_id,true,true,9);
netcdf.defVarDeflate(ncid,thickness_source_var_id,true,true,9);
netcdf.defVarDeflate(ncid,ice_var_id,true,true,9);
netcdf.defVarDeflate(ncid,mask_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,x_var_id,x);
netcdf.putVar(ncid,y_var_id,y);
netcdf.putVar(ncid,year_var_id,single(year));
netcdf.putVar(ncid,vx_var_id,ipermute(vx,[2 1]));
netcdf.putVar(ncid,vy_var_id,ipermute(vy,[2 1]));
netcdf.putVar(ncid,v_source_var_id,ipermute(v_source,[2 1]));
netcdf.putVar(ncid,thick_var_id,ipermute(thickness,[2 1]));
netcdf.putVar(ncid,thickness_source_var_id,ipermute(thickness_source,[2 1]));
netcdf.putVar(ncid,ice_var_id,ipermute(uint8(ice),[2 1 3]));
netcdf.putVar(ncid,mask_var_id,ipermute(iceshelf_mask,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp 'done'

%% 
% 
% load('icemask_composite.mat'); 
% 
% 
% filename = '/Users/cgreene/Documents/GitHub/ice-shelf-geometry/data/icemask_composite.nc'; 
% 
% disp 'saving data now...'
% 
% %iceshelf_mask = uint8(iceshelf_mask); 
% 
% % 1. Create netCDF file handle and Attributes
% mode = netcdf.getConstant('NETCDF4');
% mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
% ncid=netcdf.create(filename,mode);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Antarctic time-evolving ice mask.');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','.');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date_created',datestr(now,'yyyy-mm-dd'));
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Institution','NASA Jet Propulsion Laboratory (JPL), California Institute of Technology');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Citation','Chad A. Greene, Alex S. Gardner, Nicole-Jeanne Schlegel & Alexander D. Fraser, 2022. Antarctic calving loss rivals ice-shelf thinning. Nature. https://doi.org/10.1038/s41586-022-05037-w');
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GitHub','https://github.com/chadagreene/ice-shelf-geometry');
% 
% % 2. Define dimensions
% mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
% netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
% netcdf.putAtt(ncid,mapping_var_id,'false_northing',0);
% netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
% netcdf.putAtt(ncid,mapping_var_id,'inverse_flattening',298.257224);
% netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
% netcdf.putAtt(ncid,mapping_var_id,'scale_factor_at_projection_origin',1);
% netcdf.putAtt(ncid,mapping_var_id,'semi_major_axis',6378137);
% netcdf.putAtt(ncid,mapping_var_id,'spatial_epsg',3031);
% netcdf.putAtt(ncid,mapping_var_id,'spatial_proj','+proj=stere +lat_0=90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs');
% netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',-71);
% netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',0);
% 
% % Define x 
% x_id     = netcdf.defDim(ncid,'x',length(x));
% x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
% netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate');
% netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
% netcdf.putAtt(ncid,x_var_id,'units',        'meter');
% 
% % Define y
% y_id     = netcdf.defDim(ncid,'y',length(y));
% y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
% netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate');
% netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
% netcdf.putAtt(ncid,y_var_id,'units',        'meter');
% 
% % Define year
% year_id     = netcdf.defDim(ncid,'year',length(year));
% year_var_id = netcdf.defVar(ncid,'year','NC_FLOAT',year_id);
% netcdf.putAtt(ncid,y_var_id,'long_name',    'time');
% netcdf.putAtt(ncid,y_var_id,'units',        'decimal year');
% 
% % Define ice
% ice_var_id = netcdf.defVar(ncid,'ice','NC_BYTE',[x_id y_id year_id]);
% netcdf.putAtt(ncid, ice_var_id,'long_name','Binary ice grid');
% netcdf.putAtt(ncid, ice_var_id,'grid_mapping', 'mapping');
% 
% % Compress and stop variable definition
% netcdf.defVarDeflate(ncid,ice_var_id,true,true,9);
% netcdf.endDef(ncid);
% 
% %3. Place data
% netcdf.putVar(ncid,x_var_id,single(x));
% netcdf.putVar(ncid,y_var_id,single(y));
% netcdf.putVar(ncid,year_var_id,single(year));
% netcdf.putVar(ncid,ice_var_id,ipermute(uint8(ice),[2 1 3]));
% 
% %4. Close file 
% netcdf.close(ncid)
% 
% disp 'done'