%% main_generate_coeffs
% loops through all CTDs to make and save EOF weights
% eeshan bhatt
% 
% MIT Laboratory for Autonomous Marine Sensing Systems
% Tactical Decision Aids for ICEX20
%
% supplement to IEEE JOE submission
% DOI: -----
% 
% feel free to reach out to ebhatt (at) whoi (dot) edu if you would like to
% walk through a full demo or ask any questions about our operational
% deployment

%% prep workspace
clear; clc; close all;

% color map
mapColor = containers.Map({'RBR','XCTD'},...
    {[153 51 153]./256,[200 78 0]./256});

traceGray = [0.7 0.7 0.7 0.1];

% plot dimensions
plotDepth = [0 500];
plotWidthSSP = [1425 1465];

%% start timer
tic;

%% set source depth
zs = 30; %m

%% load ctd casts - CTD
ctd_file = './data/icex20-all-ctd-write-weights.mat';
load(ctd_file);
n = 2; %1 to 10

%% load EOF file
eeof_file = './data/eeof-itp-2013.nc';
OBJ = eb_read_eeof(eeof_file,true);

%% 1 : loop through combinations

% figure handle
figAllCombinations = figure('Name','all-combinations','renderer','painters','position',[108 108 600 800]);

% find valid depths for this CTD cast
ind_valid_depths = ~isnan(CTD(n).grid_c);
maxValidDepth = max(CTD(n).grid_z(ind_valid_depths));
adjustDepthSSP = CTD(n).grid_c;

count = 0;
for ne = 1:OBJ.num_eofs
    
    % find combinations for ne EOFs
    combos = nchoosek(1:OBJ.num_eofs,ne);
    [num_combos,~] = size(combos);
    
    % loop through combinations
    for cmb = 1:num_combos
        count = count + 1;
        eof_subset = combos(cmb,:);
        
        % find weights for valid depths
        weights = OBJ.eofs(ind_valid_depths,eof_subset) \ (adjustDepthSSP(ind_valid_depths) - OBJ.baseval(ind_valid_depths));
        new_weights = zeros(OBJ.num_eofs,1);
        new_weights(eof_subset) = weights;
        A(count).sspWeights = new_weights;
        yhat = OBJ.baseval + OBJ.eofs*new_weights;
        A(count).sspEstimate = yhat;
        A(count).sspDepth = OBJ.depth;
        
        % output nicely formatted string
        str = '-  -  -  -  -  -  -  ';
        for es = eof_subset
                str(3.*es-2) = num2str(es);
        end
        A(count).eofSubset = str;
        
        % plot ssp traces
        hold on
        plot(yhat,OBJ.depth,'color',traceGray);
    end
    plot(CTD(n).bin_c,CTD(n).bin_z,'-','color',[mapColor(CTD(n).type) .1]);
    plot(CTD(n).grid_c,CTD(n).grid_z,'.','color',mapColor(CTD(n).type),'MarkerSize',15);
    hold off
    grid on
    set(gca,'ydir','reverse');
    title(['CTD: ' CTD(n).type  ' at ' datestr(CTD(n).time,'yyyy mmm dd HHMM')],'fontsize',18,'units','normalized','position',[0.5 1.03 0]);
    xlim(plotWidthSSP);
    xlabel('c [m/s]');
    ylim(plotDepth);
    ylabel('depth [m]');
end
text(max(plotWidthSSP),-.005*range(plotDepth),['All EOF combinations N=' num2str(count)],'HorizontalAlignment','right','VerticalAlignment','bottom')

%% 2 : plot KDEs

figKDE = figure('Name','kde','renderer','painters','position',[215 690 1400 240]);

plotWeights = [A.sspWeights];

for ne = 1:OBJ.num_eofs
    subplot(1,OBJ.num_eofs,ne);
    
    xline(0,'--','color',mapColor(CTD(n).type),'linewidth',2);
    hold on
    
    % plot found weights
    thisWeight = unique(round(plotWeights(ne,:),2));
    for tw = thisWeight
        plot([tw tw],[0 1],'-','color',[mapColor(CTD(n).type) 0.1],'linewidth',1);
    end

    % plot kde
    plot(OBJ.xi(ne,:),OBJ.f(ne,:)/max(OBJ.f(ne,:)),'color','k');
    grid on
    xval = max(abs(OBJ.xi(ne,:)));    
    
    hold off
    % plot pretty
    if ne > 1
        yticklabels([]);
    else
        ylabel('pdf');
    end
    
    xlabel(sprintf('w_%d',ne))
    set(gca,'fontsize',12)
    
    xPost = max(abs([xval thisWeight]));
    xPost = xPost + 0.1.*range([-xPost xPost]);
    xlim([-xPost xPost]);
end

% overall plot things
sgtitle('Weight distribution from the training set','fontsize',14,'fontweight','bold');


%% 3 : calculate environmental error

% figure handle
figEnvError = figure('Name','env-error','renderer','painters','position',[108 108 600 800]);

for c = 1:count
    ydiff = A(c).sspEstimate - CTD(n).grid_c;
    ydiff = ydiff(ind_valid_depths);
    A(c).envDepthError = ydiff;
    A(c).envMeanError = mean(abs(ydiff));
    A(c).envStdError = std(abs(ydiff));
    A(c).envDepth   = A(c).sspDepth(ind_valid_depths);
end

% only really look into top half
meanEnvError = mean([A.envMeanError]);
if isnan(meanEnvError)
    error('environmental error returns NaN')
end

for c = 1:count
    if A(c).envMeanError <= meanEnvError
        A(c).envStatus = 1;
    else
        A(c).envStatus = 0;
    end
end
envErrorIndValid = sum([A.envStatus]);

% figure
for c = 1:count
    
    hold on
    
    if A(c).envStatus
        plot(A(c).envDepthError,A(c).envDepth,'color',[mapColor(CTD(n).type) 0.05]);
    else
        plot(A(c).envDepthError,A(c).envDepth,'color',traceGray);
    end
    
end
scatter(zeros(size(CTD(n).grid_z)),CTD(n).grid_z,20,'filled','MarkerFaceColor',mapColor(CTD(n).type));
hold off
grid on
set(gca,'ydir','reverse')
title('Depth-dependent sound speed error','fontsize',18,'units','normalized','position',[0.5 1.03 0]);
xlabel('\deltac [m/s]');
ylim([0 max(CTD(n).bin_z)]);
ylabel('CTD depth [m]');
text(max(xlim()),-.005*range(ylim()),['CTD: ' CTD(n).type  ' at ' datestr(CTD(n).time,'yyyy mmm dd HHMM')],'HorizontalAlignment','right','VerticalAlignment','bottom')

%% calculate ray shift

% pin the tail on the sound speed
max_ctd_depth = max(CTD(n).bin_z);
ind_eof_depth = find(OBJ.depth <= max_ctd_depth,1,'last');
ind_ctd_depth = find(CTD(n).bin_z == OBJ.depth(ind_eof_depth));
ssp_offset = CTD(n).bin_c(ind_ctd_depth) - OBJ.baseval(ind_eof_depth);
ssp_tail = ssp_offset.*ones(size(OBJ.baseval(ind_eof_depth+1:end))) + OBJ.baseval(ind_eof_depth+1:end);
depth_tail = OBJ.depth(ind_eof_depth+1:end);
ssp_merge = [CTD(n).bin_c(1:ind_ctd_depth) ssp_tail.'];
depth_merge = [CTD(n).bin_z(1:ind_ctd_depth) depth_tail.'];

% set up for ray tracer
zMax = 2680;
Zq = 0:1:zMax;
Cq = interp1(depth_merge,ssp_merge,Zq);

theta0 = [-30:5:-15 -10:2:10 15:5:30];
ttSpread = [0.6 1.3 2 3 4 5]; % s
sstep = 1; % m
ssp_zs = interp1(Zq,Cq,zs);
numstep = round (1.5 .* sstep .* ssp_zs .* max(ttSpread) );
[ctdRo,ctdZo,ctdTo] = eb_raytrace(zs,theta0,numstep,sstep,Zq,Cq,0,zMax);
[Ro,Zo] = h_interpolate_by_owtt(ttSpread,ctdRo,ctdZo,ctdTo);

Ri = NaN(numel(theta0),count,numel(ttSpread));
Zi = Ri;

fprintf('Raytracing:');
fprintf(['\n' repmat('.',1,envErrorIndValid) '\n\n']);
parfor c = 1:count
        
    if A(c).envStatus
        
        % fix grid spacing of SSP for ray tracer
        Cq = interp1(A(c).sspDepth,A(c).sspEstimate,Zq);
        
        % ray tracer -- sometimes EOFs break ray tracing, keep track
        try
            fprintf('\b|\n');
            [R,Z,T] = eb_raytrace(zs,theta0,numstep,sstep,Zq,Cq,0,zMax);
            A(c).raytraceR = R;
            A(c).raytraceZ = Z;
            A(c).raytraceT = T;
            % accumulate eigenray shift --- ttSpread x #EOFcombinations x thetas
            [Ri(:,c,:),Zi(:,c,:)] = h_interpolate_by_owtt(ttSpread,R,Z,T);
        catch
            A(c).envStatus = 0;
        end
    end
end

rayErrorInvalid = sum([A.envStatus]);
fprintf('caught ray errors on = %u \n',envErrorIndValid - rayErrorInvalid);

%% 4 : plot rayshift on raytrace

% figure handle
figAbsoluteRayShift = figure('Name','ray-shift-on-raytrace','renderer','painters','position',[108 108 1000 700]); clf;

% add ray trace
invalid = ctdTo > max(ttSpread);
plotR = ctdRo; plotZ = ctdZo;

plotR(invalid) = NaN;
plotZ(invalid) = NaN;
plot(plotR.',plotZ.','color',[mapColor(CTD(n).type) 0.5],'linewidth',1.5,'handlevisibility','off');
colorMax = max(ttSpread)+0.5;

hold on
for t = 1:numel(ttSpread)
    % add point cloud
    for c = 1:count
        
        if A(c).envStatus
            cloudR = squeeze(Ri(:,c,t));
            cloudZ = squeeze(Zi(:,c,t));
            
            scatter(cloudR,cloudZ,30,'filled','MarkerFaceColor',ones(3,1).*(colorMax-ttSpread(t))./colorMax,'MarkerFaceAlpha',0.4,'handlevisibility','off');
        end
    end
    
    scatter(NaN,NaN,30,'filled','MarkerFaceColor',ones(3,1).*(colorMax-ttSpread(t))./colorMax,'MarkerFaceAlpha',0.9);
    lgdStr{t} = sprintf('t = %1.1f s',ttSpread(t));
end
hold off

% make plot better
set(gca,'ydir','reverse');
grid on
ylim([0 zMax]);
xlim([0 max(plotR(:))+100]);
title('Ray shift in the water column','fontsize',18)
ylabel('depth [m]');
xlabel('range [m]');
lgd = legend(lgdStr,'location','NorthEastOutside');

%% 5 : plot rayshift agnostic to <r,z>

% figure handle
figRelativeRayShift = figure('Name','ray-shift','renderer','painters','position',[108 108 1000 700]); clf;
gvel = interp1(CTD(n).bin_z,CTD(n).bin_c,zs);

for t = 1:numel(ttSpread)
    
    subplot(2,3,t);
    
    deltaR = [];
    deltaZ = [];
    for c = 1:count
        if A(c).envStatus
            A(c).rayRangeError(t,:) = (squeeze(Ri(:,c,t)) - Ro(:,t));
            A(c).rayDepthError(t,:) = (squeeze(Zi(:,c,t)) - Zo(:,t));
            A(c).rayMeanError = mean(sqrt(A(c).rayRangeError(t,:).^2 + A(c).rayDepthError(t,:).^2));
            A(c).rayStdError = std(sqrt(A(c).rayRangeError(t,:).^2 + A(c).rayDepthError(t,:).^2));
            
            % make this relative for the plot
            deltaR = [deltaR; (squeeze(Ri(:,c,t)) - Ro(:,t))];
            deltaZ = [deltaZ; (squeeze(Zi(:,c,t)) - Zo(:,t))];
        end
    end
    
    % make these usefully relative
    deltaR = deltaR ./ (gvel .* ttSpread(t))  .*100;
    deltaZ = deltaZ ./ zMax .* 100;
    
    % get colors for scatter based on how many things are "close" to it
    alphaVal = zeros(1,numel(deltaR));
    for d = 1:numel(deltaR)
        indR = abs(deltaR - deltaR(d) > 2.*std(deltaR));
        indZ = abs(deltaZ - deltaZ(d) > 2.*std(deltaZ));
        alphaVal(d) = sum(indR) + sum(indZ);
    end
    alphaVal(alphaVal <= mean(alphaVal) - std(alphaVal)) = mean(alphaVal) - 0.5*std(alphaVal); 
    alphaVal = alphaVal/max(alphaVal) -.1;
    alphaVal = alphaVal .';
    
    hold on
    scatter(deltaR,deltaZ,30,alphaVal.*ones(1,3),'filled','MarkerFaceAlpha',traceGray(4));
    scatter(0,0,60,'filled','MarkerFaceColor',mapColor(CTD(n).type));
    hold off
    
    % things for all plots
    title(sprintf('t = %1.1f s',ttSpread(t)),'fontweight','normal','fontsize',14);
    grid on
    
    %axis equal
    set(gca,'fontsize',12)
    
    if mod(t,3) == 1
        ylabel('% depth shift');
    end
    
    if t >= 4
        xlabel('% range shift');
    end
end
sgtitle('Relative ray shift','fontsize',18,'fontweight','bold')

%% end timer 
toc;

%% make table for interactivity!

% get index
simValid = find([A.envStatus]==1).';
% eof_subset
stat.subset = {A(simValid).eofSubset}.';
% ssp error
stat.sspError = [A(simValid).envMeanError].';
% ssp std
stat.sspSTD = [A(simValid).envStdError].';
% ray error
stat.rayError = [A(simValid).rayMeanError].';
% ray std
stat.raySTD   = [A(simValid).rayStdError].';
% rank
errorPenalty  = zscore(stat.sspError) + zscore(stat.rayError) + zscore(stat.raySTD);
[errorPenalty,orderEP]   = sort(errorPenalty);

% re-order
stat.subset = stat.subset(orderEP);
stat.sspError = stat.sspError(orderEP);
stat.rayError = stat.rayError(orderEP);
stat.raySTD   = stat.raySTD(orderEP);
stat.rank     = [1:numel(orderEP)].';
stat.indexEOF = simValid(orderEP);

% make table
interactiveTable = struct2table(stat);
interactiveTable = movevars(interactiveTable,'rank','before','subset');

% write table
writetable(interactiveTable,sprintf('ctd%u-table.txt',n),'Delimiter',',');

% display table
display(interactiveTable);

%% loop through interactivity
% chosen weight gets highlighted in blue

% myBlue
myBlue = [51, 152, 152]./256;

% while loop for choosing
chooseRank = -1;
countRank = 0;
fprintf('\nPlease enter a rank to see it projected on all plots [enter 0 to save weights for your last rank] \n')
while chooseRank~=0
    
    countRank = countRank + 1;
    
    chooseRank = input('      rank: ');
    if chooseRank >= 1
        
        % match back to EOF realization case
        iso = stat.indexEOF(chooseRank);
        
        %% ----- figAllCombinations
        figure(figAllCombinations);
        
        % manage text on this plot only
        if countRank == 1
            t0charArray = char('\it 1','\it 2','\it 3','\it 4','\it 5','\it 6','\it 7');
            t0 = text(1459,20,t0charArray,'VerticalAlignment','top','HorizontalAlignment','left','color',myBlue);
        else
            delete(t1);
            delete(t2);
        end
        
        % remove last rank choice projection
        h_clean_plot(countRank,1);
        
        % plot sspEstimate
        hold on
        h1 = plot(A(iso).sspEstimate,A(iso).sspDepth,'color',myBlue);
        hold off
        t1 = text(max(plotWidthSSP),20,num2str(A(iso).sspWeights,'%2.3f'),'VerticalAlignment','top','HorizontalAlignment','right','color',myBlue);
        t2 = text(1459,18,sprintf('rank = %u',chooseRank),'VerticalAlignment','bottom','HorizontalAlignment','left','fontweight','bold','fontsize',12,'color',myBlue);
        
        %% ----- figKDE
        figure(figKDE);
        chooseWeights = A(iso).sspWeights;
        
        for ne = 1:OBJ.num_eofs
            subplot(1,OBJ.num_eofs,ne);
            h_clean_plot(countRank,1);
            
            hold on
            plot([chooseWeights(ne) chooseWeights(ne)],[0 1],...
                'color',myBlue,'linewidth',2)
            hold off
        end
        %% ----- figEnvError
        figure(figEnvError);
        h_clean_plot(countRank,1);
        
        hold on
        plot(A(iso).envDepthError,A(iso).envDepth,'-o','color',myBlue);
        hold off
        
        %% ----- figAbsoluteRayShift
        figure(figAbsoluteRayShift)
        
        if countRank > 1
            delete(p3);
        end
            
        miniRayR = A(iso).raytraceR;
        miniRayZ = A(iso).raytraceZ;
        miniRayT = A(iso).raytraceT;
        
        tooLong = miniRayT > max(ttSpread);
        miniRayR(tooLong) = NaN;
        miniRayZ(tooLong) = NaN;
        
        hold on
        p3 = plot(miniRayR.',miniRayZ.','color',[myBlue 0.3],'handlevisibility','off');
        hold off
        
        %% ----- figRelativeRayShift
        figure(figRelativeRayShift)
        
        for t = 1:numel(ttSpread)
            subplot(2,3,t)
            
            h_clean_plot(countRank,2);
            
            miniDeltaR = A(iso).rayRangeError(t,:) ./ (gvel .* ttSpread(t)) .* 100;
            miniDeltaZ = A(iso).rayDepthError(t,:) ./ zMax .* 100;
            
            hold on
            scatter(miniDeltaR,miniDeltaZ,30,'filled','MarkerFaceColor',myBlue,'MarkerFaceAlpha',0.7);
            scatter(0,0,30,'filled','MarkerFaceColor',mapColor(CTD(n).type));
            hold off
        end
        
    coeffs.chooseRank = chooseRank;
    %% export weights
    elseif chooseRank == 0
        fprintf('\nthank you : chosen weights ... \n\n');
        for k = 1:numel(A(iso).sspWeights)
            fprintf('w%u = %3.3f \n',k,A(iso).sspWeights(k));
        end
        
        % save all figures
        figure(figAllCombinations);
        h_printPNG(sprintf('ctd%u-allCombinations.png',n));
        
        figure(figKDE);
        h_printPNG(sprintf('ctd%u-KDE.png',n));
        
        figure(figEnvError);
        h_printPNG(sprintf('ctd%u-envError.png',n));
        
        figure(figAbsoluteRayShift);
        h_printPNG(sprintf('ctd%u-absoluteRayShift.png',n));
        
        figure(figRelativeRayShift);
        h_printPNG(sprintf('ctd%u-relativeRayShift.png',n));
        
        % write to mat file
        coeffs.eofWeights = A(iso).sspWeights;
        coeffs.timeWritten = datestr(now());
        save(sprintf('ctd%u-coeffs.mat',n),'coeffs');
    end
end

%% clean up
close all; clear;

%% helper function : h_interpolate_by_owtt
function [Ri,Zi] = h_interpolate_by_owtt(ttSpread,R,Z,T)

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

%% helper function: h_clean_plot
function [] = h_clean_plot(count,num2remove)
if count >= 2
    children = get(gca,'children');
    delete(children(1:num2remove));
end
end

%% helper function: h_printPNG
function [filename] = h_printPNG(filename)
% export_fig(filename,'-q101','-r300','-painters','-png');

% saveas(gcf,[filename '.svg']);
end