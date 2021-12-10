%% Helper Function: getCVal;

function [c,cz,czz]=h_getCVal(r,z,cVal)
%return a c value in m/s for a given r,z value:
%search grid
c=0; % soundspeed
cz=0; % slope
czz=0; % second derivative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TO GET soundspeed and soundspeed gradient for r,z  %%%%%%%%%%%%

% cVal is range-independent, don't need r;

% if ~isreal(z)
%     warning('Depth is imaginary:')
%     display(z)
% end

c = cVal{1}(z);
cz = cVal{2}(z);
czz = cVal{3}(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end