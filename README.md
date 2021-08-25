# ice-shelf-geometry
Matlab scripts that use all available observations to grow and trim the extents of Antarctic ice shelves.

## Script Workflow 
1. **`iceshelf_mask_generator.m`** uses Mouginot's iceshelves\_2008\_v2 outlines to create `iceshelf_mask.mat`, which contains a 240 m resolution mask on the ITS\_LIVE. This script also dilates the ice shelf mask by 100 km to account for any possible ice shelf growth. 
2. **`flow_dem_extend.m`** creates reference ice thickness and velocity grids that extend 100 km beyond the ~2008 extents of the ice sheet. Velocity is obtained by combining ITS\_LIVE and MEaSUREs v2 data as an error-weighted average of the two datasets, then velocities are extrapolated as constant values ~100 km beyond the calving front, in the direction of ice flow. Ice thickness is inverted from BedMachine v2 surface elevations (which are REMA - FAC). Where calving occurred before the ~2014? nominal date of BedMachine and BedMachine contains 0 thickness, surface elevations are taken from REMA, then Bedmap2, then Bamber DEM, then RAMP2, and from each of these non-BedMachine surface elevations, FAC from GEMB is subtracted.  

## Functions and other files 
* **`inpolygon_map()`** is a _much_ faster version of `inpolygon` used by `iceshelf_mask_generator.m` to mask ice shelf boundaries. 
* **iceshelves\_2008\_v2\_names.csv** contains the names and numerical indices of all 181 ice shelves in the extended ice shelf mask. 

## Datasets created by these scripts
The datasets created by these scripts are too large to be included in this GitHub repository.  

