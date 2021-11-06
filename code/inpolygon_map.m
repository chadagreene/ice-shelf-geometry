function in = inpolygon_map(x,y,xi,yi) 
% inpolygon_map is a much faster implementation of inpolygon, for strictly
% gridded coordinates that are equally spaced. 
% 
%% Syntax 
% 
%  in = inpolygon_map(x,y,xi,yi)
% 
%% Description 
% 
% in = inpolygon_map(x,y,xi,yi) creates an MxN logical mask "in" of a grid
% whose spatial coordinates are given by x,y. If x,y are 1D vectors,
% the length of y corresponds to the number of rows in the output mask, and
% the length of x corresponds to the number of columns in the output mask.
% Inputs xi,yi specify the boundaries of any polygons, and can be arrays or
% NaN-separated arrays. If xi,yi are NaN-separated arrays, counterclockwise
% polygons will be true, and clockwise polygons will be holes. 
% 
%% Requirements 
% Requires Image Processing Toolbox. And if you want to enter nan-separated
% arrays xi,yi, then you'll also need the Mapping Toolbox. 
%
%% Author Info
% Written by Chad A. Greene of NASA Jet Propulsion Laboratory, 
% August 2021. 
% 
% See also: inpolygon, poly2mask, and roipoly. 

%% Input checks 

assert(license('test','image_toolbox'),'Sorry, the inpolygon_map function requires the image processing toolbox.') 

narginchk(4,4) 
if ~isvector(x) 
   assert(isequal(size(x),size(y)),'Inputs x and y can only be equally-spaced vectors, or grids as if created by meshgrid.') 
   x = x(1,:); % assumes x and y were created by [x_grid,y_grid] = meshgrid(x_vector,y_vector); 
   y = y(:,1); 
end
assert(isequal(size(xi),size(yi)),'Error: Inputs xi,yi must have matching dimensions.') 

dx = diff(x); 
dy = diff(y); 
assert(all(dx==dx(1)),'Error: input gridpoints x,y must be equally spaced.')
assert(all(dy==dy(1)),'Error: input gridpoints x,y must be equally spaced.')

%%

x_extent = [x(1) x(end)]; 
y_extent = [y(1) y(end)]; 

% Transform xi,yi into pixel coordinates.
roix = axes2pix(numel(x), x_extent, xi);
roiy = axes2pix(numel(y), y_extent, yi);

if any(isnan(roix))
   assert(license('test','image_toolbox'),'Sorry, if xi,yi contain any NaNs, the inpolygon_map function requires the mapping toolbox.') 

   [roixc,roiyc] = polysplit(roix,roiy); % removes nans 
   isccw = ~ispolycw(roixc,roiyc); % clockwise usually corresponds to holes. 

   in = false(numel(y),numel(x)); % preallocates 
   holes = in; 

   % Loop through each section: 
   for k = 1:length(roixc)
      if isccw(k)
         in(poly2mask(roixc{k}, roiyc{k}, numel(y), numel(x))) = true; 
      else
         holes(poly2mask(roixc{k}, roiyc{k}, numel(y), numel(x))) = true;
      end
   end
   
   % Remove holes: 
   in(holes) = false;  

else
   in = poly2mask(roix, roiy, numel(y), numel(x));
end

end
