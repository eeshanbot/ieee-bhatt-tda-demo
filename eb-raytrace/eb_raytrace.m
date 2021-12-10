function [Ri,Zi,Ti] = eb_raytrace(zsource,theta,numstep,sstep,sspDepthGrid,sspVal,bathyRangeGrid,bathyVal)
% modified from 2.681 Env. Ocean Acoustics HW 3
%
% zsource           the source depths you want for beams
% theta             the set of beam angles you want to run the raytrace code at,
%                       between -90 and 90 degrees
% numstep           the number of steps to run
% sstep             the step size along the ray
% sspDepthGrid      the depth grid of the sound speed profile
% sspVal            the sound speed profile
% bathyRangeGrid    the range grid of the bathymetry [if range independent, put 0]
% bathyVal          the bathymetry [if range indepenent, put one value]

%% Interpolate Data
% If you've read COA, you'll know that you will need spatial derivatives of
% the sound speed field. It behooves you to calculate those derivatives
% first, and access them as nececessary, instead of making an interpolation
% every time.

theta = theta * -pi/180;

[cVal,dVal] = h_dataInterp(sspDepthGrid,sspVal,bathyRangeGrid,bathyVal);

% ray tracing parameters
A0 = 1;
z0 = zsource;

% loop over all thetas
for dth=1:length(theta) % get theta value
    %% Initialize variables needed from raytracing
    
    % get sound speed, 1st derivative of sound speed, 2nd derivative of sound speed
    [c0,cz0,czz0]=h_getCVal(0,z0,cVal);
    
    % for each theta value, calculate ray
    theta0 = theta(dth);
    
    % initialize raytrace vectors
    sray=[0];%s-values for ray
    xiray=[cos(theta0)/c0];% xi values for ray
    zetaray=[sin(theta0)/c0];% zeta values for ray
    thetavec=[theta0]; % theta values for ray
    rray=[0];%range vector for ray
    zray=[z0];%depth vector for ray
    cray=[c0];%soundspeed vector for ray
    tau=[0];%tau vector for ray
    q=[0];
    p=[1/c0];
    A=[A0]; %ray amplitude
    cz=[cz0]; %first derivative of ssp
    czz=[czz0]; %second derivative of ssp
    J0=0;
    i=1; %step count
    rsign=1;
    
    %% Iterate over numsteps to calculate out the ray
    
    while i<numstep
        
        r=rray(i);%r from previous step
        z=zray(i);%z from previous step
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% INSERT CODE HERE TO GET r, z OF RAY IN NEW STEP %%%%%%
        
        r=r+sstep*cos(thetavec(i));
        z=z+sstep*sin(thetavec(i));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %figure out if there is bottom or surface interaction,
        %calculate new values
        
        [interact,sray,zray,rray,zetaray,xiray,cray,tau,thetavec,c0,theta0,cz,czz]=h_getInteraction(r,z,...
            dVal,sray,zray,rray,zetaray,xiray,cray,theta0,cVal,tau,thetavec,c0,cz,czz);
        
        %If there is an interaction, add two to 'i' and get new p, q, A
        if interact
            
            %add 2 to i to account for the two points from interaction fctn
            [q,p,A]=h_getAmplitudes(cray,zray,rray,sray(i)-sray(i-1),xiray,theta0,A0,i,p,q,A,c0,cz,czz); %First addition
            [q,p,A]=h_getAmplitudes(cray,zray,rray,sray(i+1)-sray(i),xiray,theta0,A0,i+1,p,q,A,c0,cz,czz); %Second
            i=i+2; %get off of the surface
            
            continue
        end
        
        % add the new rray, zray values
        rray(i+1)=r; % the values you just calculated
        zray(i+1)=z;
        
        % Get new soundspeed for new location
        [cray(i+1), cz(i+1), czz(i+1)]=h_getCVal(rray(i+1),zray(i+1),cVal);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% PUT YOUR CODE TO GET zetaray(i+1) HERE %%%%%%%%%%%%%%
        
        zetaray(i+1) = zetaray(i) + -1/cray(i+1)^2 * cz(i+1)*sstep;
        xiray(i+1) = xiray(i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% PUT YOUR CODE TO GET NEW thetavec(i+1) HERE  %%%%%%%%%%
        
        thetavec(i+1) = atan2(zetaray(i+1),xiray(i+1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Move Forward
        sray(i+1)=sray(i)+sstep;
        dtau=sstep/cray(i);
        tau(i+1)=tau(i)+dtau;
        
        % next, calculate the amplitudes,p and q:
        [q,p,A]=h_getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz,czz);
        
        i=i+1;
    end
    
    % just in case the ray ends on a bounce!
    rray = rray(1:numstep);
    zray = zray(1:numstep);
    tau = tau(1:numstep);
    
    Ri(dth,:) = rray;
    Zi(dth,:) = zray;
    Ti(dth,:) = tau;
end

end
