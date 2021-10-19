These files can be used to set up ISSM experiments for to compare the effects of calving vs ice shelf thinning on grounding line flux, as measured by the instantaneous velocity response. 

All data in these files are on the ITS_LIVE 240 m grid. 

The four experiments: 
1. Observed calving 1997 to 2021: 24 ~annual masks are in icemask_composite.mat. 
2. Hypothetical calving 0 to 100% of all ice shelves: 101 masks in 1% increments are in icemask_buffered.mat. 
3. Observed thinning 1992 to 2018: Annual thickness cube is in iceshelf_thickness_cube_1992-2018.mat. 
4. Hypothetical future thinning if 1992-2018 rates continue: 50 years of future thinning from 2018, in 2 year increments in iceshelf_thickness_cube_future.mat. 

Notes: 
* In some of the files below I've named the variable `thickness`, but elsewhere it's `H`. It's the same variable, sorry for switching it up. 
* The uint16 data uses a fill value of 32767 to indicate ocean. 

FILES

issm_calving_melt_setup.mat 
* Use this to set up the ISSM model. 
* I used rho_ice=917 and rho_seawater=1027 for hydrostatic calculations. 
* NaN-out everything that is not ever_ice. 
* created by issm_calving_melt_setup.m
* Variables: 
- bed from 'BedMachineAntarctica-2021-03-03.nc'
- ever_ice logical mask showing true for any grid cell that has ever been observed to be ice (or rock, but we're calling that ice). 
- ground anything that is not ocean or ice shelf in the BedMachine mask. 
- vx,vy extended and trimmed to the ever_ice mask. 
- thickness composite extended and trimmed to ever_ice mask. 

icemask_composite.mat
* An "ice cube" of 24 ~annual ice masks, compiled MOA, RAMP, Sentinel-1, and MODIS.
* created by ice_compiler.m. 
* Variables: 
- ice 18392x22896x24 logical cube that is true wherever ice or rock was present. 
- cx,cy 1x24 cell arrays containing the contoured outline of the ice cube. 

icemask_buffered.mat 
* Ice masks buffered inland from from 0 to 100% of each ice shelf's area, buffered from maximum observed extents. 
* created by icemask_buffer.m. 
* Variables: 
- ice 18392x22896x101 a logical "ice cube" 

iceshelf_thickness_cube_1992-2018.mat
* Ice thickness cube extrapolated (and in some places trimmed) from Fernando Paolo's FULL_CUBE_v4.h5.
* Ice thickness was tied to 2014, then anomalies from Fernando's H_filt10 thickness data, and thinning that has occurred since 1992 was added back to get the 1992 thickness. 
* Thickness was limited to prevent changes in the grounding line location, using hydrostatic approximation with rho_ice=917 and rho_seawater=1027.
* Although in real life the area of the ice sheet changed due to calving and regrowth, the areal extents of this thickness cube is the same for all years. The extents are defined as the grid cells that were always ice from the years 1997-2018. 
* Fill Value 32767.
* created by iceshelf_thickness_cube_generator.m.
* Variables: 
- H ice thickness in meters observed 1992-2018. 
- always_ice is true for all grid cells that were ice every year from 1997 (the earliest calving front observation) to 2018 (the latest thickness observation). 
- ground is true for all grid cells that are not ocean or ice shelf in 'BedMachineAntarctica-2021-03-03.nc'.

iceshelf_thickness_cube_future.mat
* Thickness trends starting from the last observed thickness (2018) extrapolated into the future 50 years. 
* Thickness was limited to prevent changes in the grounding line location, using hydrostatic approximation with rho_ice=917 and rho_seawater=1027.
* created by iceshelf_thickness_cube_generator.m.
* Variables
- H thickness for the next 50 years (in two year increments) if recent trends continue.
- always_ice is true for all grid cells that were ice every year from 1997 (the earliest calving front observation) to 2018 (the latest thickness observation).
- ground is true for all grid cells that are not ocean or ice shelf in 'BedMachineAntarctica-2021-03-03.nc'.
- future_years number of years into the future, starting at the last observed thickness. 

Everything here was written by Chad Greene, NASA Jet Propulsion Laboratory, October 2021. 
