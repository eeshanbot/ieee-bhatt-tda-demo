%% helper function : eb_ray_interpolate_owtt
function [Ri,Zi] = eb_ray_interpolate_owtt(ttSpread,R,Z,T)

% number of thetas to aggregrate over
numTheta = size(R,1);
for n = 1:numTheta
    
    interp_range = interp1(T(n,:),R(n,:),ttSpread);
    interp_depth = interp1(T(n,:),Z(n,:),ttSpread);
    
    % matrix to store eigenrays -- [theta x ttSpread]
    Ri(n,:) = interp_range;
    Zi(n,:) = interp_depth;
end
end
