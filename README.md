# ice-shelf-geometry
Matlab scripts that use all available observations to grow and trim the extents of Antarctic ice shelves.

# Script Workflow 
1. `iceshelf_mask_generator.m` uses Mouginot's iceshelves\_2008\_v2 outlines to create `iceshelf_mask.mat`, which contains a 240 m resolution mask on the ITS\_LIVE. This script also dilates the ice shelf mask by 100 km to account for any possible ice shelf growth. 

# Included functions 
* `inpolygon_map` is a _much_ faster version of `inpolygon` used by `iceshelf_mask_generator.m` to mask ice shelf boundaries. 
