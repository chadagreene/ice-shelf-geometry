
% Chad Greene 
% November 2021. 
% NASA Jet Propulsion Laboratory. 

load('icemask_composite.mat','cx','cy','year')

%%

header1 = '# This is an Antarctic coastline corresponding to the decimal year in the filename.';
header2 = '# The columns are EPSG:3031 (ps71) eastings and northings in meters.';
header3 = '# These coastlines correspond to the 240 m ice mask in icemask_composte.mat.';
header4 = '# Regions of ice are separated by a row containing NaN NaN.';
header5 = '# Created by Chad A. Greene, NASA Jet Propulsion Laboratory, November 2021.'; 

for k = 1:24
   
   filename = ['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/coastline_data/antarctic_coastline_',num2str(year(k)),'.txt']; 

   writematrix(header1,filename)
   writematrix(header2,filename,'writemode','append')
   writematrix(header3,filename,'writemode','append')
   writematrix(header4,filename,'writemode','append')
   writematrix(header5,filename,'writemode','append')
   writematrix([round(cell2nancat(cx(k))) round(cell2nancat(cy(k)))],filename,'writemode','append')
end

