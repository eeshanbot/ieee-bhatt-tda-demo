%% Helper Function: getAmplitudes

function [q,p,A]=h_getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz,czz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% getAmplitudes function: get q, p, A values for ray input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize:
q(i+1)=0;
p(i+1)=0;
A(i+1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INSERT YOUR CODE HERE TO GET p, q, A values   %%%%%%%%%%%%%%%%%

% Quickly assign p(i+1) = p(i) so we can solve q(i+1) using "current"
% indices. 
q(i+1) = q(i) + cray(i)*p(i)*sstep;

% Now, reevaluate p(i+1);
% dp/ds = -cnn/c^2(s) * q(s)
cnn = cray(i+1).^2 .* czz(i+1) .* xiray(i+1).^2;
p(i+1) = p(i) - cnn ./ cray(i+1)^2 * q(i+1) * sstep;

% Now find A(i+1)
A(i+1) = A0/(4 .* pi) .* ...
    sqrt(abs(cray(i+1)*cos(theta0)./(c0 .* rray(i+1) .* q(i+1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end