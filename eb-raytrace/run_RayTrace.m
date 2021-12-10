% EeShan Bhatt
% eesh@mit.edu
% 2.681 Env. Ocean Acoustics HW 4

% Set up workspace 
close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%
% The form of the function for running the ray tracing is:
% raytrace(sstep,depmat,cmat,thetas,zsources,f0,numsteps,debug)
% sstep: the step size along the ray
% depmat: the depth matrix indexed by range, with form [r1,dep1;r2,dep2...]
% cmat: the soundspeed profile, with form [z1 c1;z2 c2...]
% thetas: the set of beam angles you want to run the raytrace code at,
% between -90 and 90 degrees
% zsources: the source depths you want for beams
% f0: frequency
% numsteps: is the number of steps to run
% debug is the flag that, if true, creates plots with each top/bottom interaction

%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you want a different soundspeed profile, the format is [z1,c1;z2,c2;...]
% for example, c_munk is constructed as:
% depvec=(1:5000)';
% zbar = 2.*(z-1300)./1300;
% cvec=1500.*(1+.00737.*(zbar-1+exp(-zbar)));
% cmunk = [depvec,cvec];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings-
% size of steps along ray in m
sstep=50;
% frequency
f0 = 100;
load c_munk

%% Part 1A: For a source depth 500m, plot rays and intensity along the beams?
load depmat_flat; %flat bottom, contains depmat
numsteps=1500;

% 1a: For a source at 500m depth, thetas=[30,20,10,5,0,-5,-10,-20,-30]
thetas=[30,20,10,5,-5,-10,-20,-30]*pi/180;
zsources = [500];
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *Describe the behavior of these rays going out to at least 5 km*

%%
% 
% These rays propogate in a Munk profile. Rays with steep departure angles
% experience bottom and top interactions. Ones with shallower departure
% angles propogate through the medium without any reflections, and turn at
% varying depths. This is because these rays abide by Snell's law such that
% c0*cos(theta0) = cf*cos(thetaf) and the changing theta0 will account for
% varying turning sound speeds.
%  
% 


%% Part 1B: For a source depth at 500m depth, what range of initial beam angles do not reflect to either the surface or bottom?

numsteps=1500;
thetas=[13 12 5 -5 -12 -13]*pi/180;
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *For a source at 500 m depth, what range of initial beam angles do not
% reflect either the surface or the bottom?*
 
%%
% This plot shows thetas of +/- 13, 12, and 5. At 13 degrees there is
% bottom and top interaction; any smaller, or +12 to -12, and there are no reflected or
% interactions. This angle is purely derived from ray tracing and only
% deals with the sound speed profile! This can be found numerically through
% Snell's Law, theta = acos(spd@source/spd@surface) = 12.8258 degrees.

%% Part 1C: What range of initial beam angles pass through a receiver 20 km away at 500 m depth?
 
numsteps = 500;
thetas=[-12:4:12]*pi/180;
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *For a source at 500m depth, what range of initial beam angles pass
% through or near to a receiver at 500m depth 20 km away?*

%% 
% A receiver in this location is in a dark shadow zone if we ignore any
% rays that have bottom or top interactions. This is a bit of a
% simplification because rays that bounce once or twice will still have
% enough intensity to be picked up by a receiver. In more robust ray
% tracing codes, such as Bellhop, the a standard value for the number of
% bounces is 4 or 5. This is because after this many bounces the rays are
% too computationally intensive for the amount of energy they provide.
%

%% Part 1D: How does the intensity of these beams change with source depth?

% 1d: Plot beam intensity for thetas=[5] and zsources=[250,500,750,1000]
numsteps = 2000;
thetas=[5]*pi/180;
zsources = [250,500,750,1000];
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *How does the intensity of these beams change with source depth?*

%%
% 
% The intensity of these beams versus source depth stays fairly constant
% until 40 km. The spikes in the intensity refer to convergence zones in
% the ray path. These convergence zones are different for different
% sources. 
% 



%% Part 2A: For a source at 500m depth, plot rays and intensity along the rays.
load depmat_slant; %slanted bottom, contains depmat
numsteps=1500;

% 2a: For a source at 500m depth, thetas=[30,20,10,5,0,-5,-10,-20,-30]
thetas=[45,30,20,10,5,0,-5,-10,-20]*pi/180;
zsources = [500];
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *Describe any interesting behavior your observe*
%%
% 
% There is lots of interesting behavior in this plot!! The 
% ranges of incident angles that do not have any bottom or top
% interactions. The most interesting part is how some rays react to
% encountering the slant. The steepest ray bounces over the slant
% but the second steepest actually propagates backward. It would be
% interesting to see how one could use this information to model seamounts
% and their elevation and location. 

% 2b: Create a zoomed in figure of the bottom interaction of one of the
% rays!

%% Part 2C: What range of initial beam angles do not reflect to either the surface or the bottom?

% 2b: Vary parameters above to answer question
thetas=[5 4 0 -4 -5]*pi/180;
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)

%%
% *For a source at 500m depth, what range of initial beam angles do not 
% reflect to either the surface or the bottom?*

%%
% Using Snell's Law (as described above in 1B/C), one can see that the initial
% angle to turn around before 3000m is 4.8675 degrees. This plot shows initial
% angles of 4 and 5 degrees - one can see that the 5 degree initial angle
% just touches the seamount whereas the 4 degree does not.
% 
% 

%% Part 3: Caustic
thetas=[-20:10:20]*pi/180;
load depmat_flat;
raytrace(sstep,depmat,cmunk,thetas,zsources,f0,numsteps,false)
