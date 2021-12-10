%% Helper function: dataInterp

function [cVal,dVal] = h_dataInterp(sspDepthGrid,sspVal,bathyRangeGrid,bathyVal)

% create sound speed interpolation objects to avoid redundant calculations
depth = sspDepthGrid;
ssp = sspVal;
cz_space = gradient(ssp, depth);
czz_space = gradient(cz_space, depth);

Fc = griddedInterpolant(depth, ssp, 'spline', 'spline');
Fcz = griddedInterpolant(depth, cz_space, 'spline', 'spline');
Fczz = griddedInterpolant(depth, czz_space, 'spline', 'spline');

cVal = {Fc Fcz Fczz};

% create depth interpolation object to avoid redundant calculations
range = bathyRangeGrid; 
bathy = bathyVal;

if numel(range) == 1
    range = [range range+100000000];
    bathy = [bathy bathy];
end

bathyr = gradient(bathy, range);

Fd = griddedInterpolant(range, bathy, 'spline', 'spline');
Fdr = griddedInterpolant(range, bathyr, 'spline', 'spline');
dVal = {Fd Fdr};
end