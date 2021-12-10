%% Helper Function: getDepth

function [D,ang]=h_getDepth(r,dVal)

D=-100;ang=0;
% r- range value in m
% dVal- cell array of depth griddedInterpolants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TO GET DEPTH D and bottom angle ang at cur r, dVal  %%%%%%%%%%

D = dVal{1}(r);
ang = atan(dVal{2}(r));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end