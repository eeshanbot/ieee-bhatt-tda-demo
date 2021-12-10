load c_munk; % contains cmunk
load depmat_flat; %slanted bottom, contains depmat
numstep=10000;

% For a source at 500m depth ... 
sstep = 1;
theta=[0];
zsource = [500];

% with bathymetry
%[R,Z,T] = eb_raytrace(zsource,theta,numstep,sstep,cmunk(:,1),cmunk(:,2),depmat(:,1),depmat(:,2));

% with flat bottom
[R,Z,T] = eb_raytrace(zsource,theta,numstep,sstep,cmunk(:,1),cmunk(:,2),0,3000);

% figure
figure(1);
plot(R.',Z.');
set(gca,'ydir','reverse');
ylim([0 max(cmunk(:,1))]);
