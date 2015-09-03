

%% REQ #1
% Video example
% Subplot(121),rearing as a skeleton 
% subplot(122),trajectory of rears through the phase space of
% ethologically relavant features

Trial = MTATrial('jg05-20120317');

%% Initialize Features
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');
fvel = xyz.vel(1,[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0) = .1;
fvel.data = log10(fvel.data);

ang = create(MTADang,Trial,xyz);
fang = create(MTADang,Trial,fxyz);

load(fullfile(Trial.spath,...
    [Trial.filebase '-walk_fet_ppc.mat']))
man = Trial.xyz.copy;
man.data = mag;
man.filter('ButFilter',3,1.5,'low');

fet = Trial.xyz.copy;
fet.data = [fvel.data,man.data];

%% Open Video File

Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-walk_ex001.avi']),'Uncompressed AVI');
Rex.open;

%% Set up Figure
hfig = figure;
hax = axes;

%% Plot JPDF of speed vs traj vect dir ppc
aind = Trial.stc{'a'};
hist2(fet(aind,:),linspace(-.8,2,100),linspace(-.2,1,100))
caxis([0,200]);
title('Trajectory of Walk Periods');
xlabel('log10 speed_{BL} (cm/s)');
ylabel('ppc_{bodydir} (AU)');



ind = Trial.stc{'w'}(s+1,:);
rtraj = {};
rstart = {};
rstop = {};
s = 1;
    rtraj{s} = animatedline;
    rtraj{s}.Parent = hax;
    rtraj{s}.Marker = '.';
    rtraj{s}.MarkerEdgeColor = 'm';
    rtraj{s}.MarkerFaceColor = 'm';
    rtraj{s}.MarkerSize = 12;
    rtraj{s}.addpoints(fet(ind(1),1),fet(ind(1),2));

    rstart{s} = line(fet(ind(1),1),fet(ind(1),2));
    rstart{s}.Parent = hax;
    rstart{s}.Marker = '^';
    rstart{s}.MarkerFaceColor = 'g';
    rstart{s}.MarkerEdgeColor = 'g';
    rstart{s}.MarkerSize = 12;
    rstart{s}.Visible = 'on';

    rstop{s} = line(fet(ind(2),1),fet(ind(2),2));
    rstop{s}.Parent = hax;
    rstop{s}.Marker = 'v';
    rstop{s}.MarkerFaceColor = 'r';
    rstop{s}.MarkerEdgeColor = 'r';
    rstop{s}.MarkerSize = 12;
    rstop{s}.Visible = 'on';
    
legend({'walking trajectory','walking onset','walking offset'},'location','NorthWest');
delete(rtraj{s});
delete(rstart{s});
delete(rstop{s});

for s = 1:400,
    rind = Trial.stc{'w'}(s+1,:);
    ind = rind + [-20,20];
    ind = ind(1):10:ind(2);

    rtraj{s} = animatedline;
    rtraj{s}.Parent = hax;
    rtraj{s}.Marker = '.';
    rtraj{s}.MarkerEdgeColor = 'm';
    rtraj{s}.MarkerFaceColor = 'm';
    rtraj{s}.MarkerSize = 12;
    rtraj{s}.addpoints(fet(ind(1),1),fet(ind(1),2));

    rstart{s} = line(fet(rind(1),1),fet(rind(1),2));
    rstart{s}.Parent = hax;
    rstart{s}.Marker = '^';
    rstart{s}.MarkerFaceColor = 'g';
    rstart{s}.MarkerEdgeColor = 'g';
    rstart{s}.MarkerSize = 12;
    rstart{s}.Visible = 'off';

    rstop{s} = line(fet(rind(2),1),fet(rind(2),2));
    rstop{s}.Parent = hax;
    rstop{s}.Marker = 'v';
    rstop{s}.MarkerFaceColor = 'r';
    rstop{s}.MarkerEdgeColor = 'r';
    rstop{s}.MarkerSize = 12;
    rstop{s}.Visible = 'off';


    for i = ind,
        rtraj{s}.addpoints(fet([i,i+10],1),fet([i,i+10],2));
        if i >= rind(1), rstart{s}.Visible = 'on'; end
        if i >= rind(2), rstop{s}.Visible = 'on'; end
        drawnow
        frm = getframe(hfig);
        Rex.writeVideo(frm.cdata);
    end

    rstart{s}.MarkerSize = 6;
    
    rtraj{s}.MarkerSize = 6;
    rtraj{s}.MarkerEdgeColor = [.6,.6,.6];
    rtraj{s}.MarkerFaceColor = [.6,.6,.6];

    rstop{s}.MarkerSize = 6;
    if s >10,
        delete(rstart{s-9})
        delete(rtraj{s-9})
        delete(rstop{s-9})
    end

end


Rex.close;

%% REQ #2


xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5,'low');
nfet = MTADfet(Trial.spath,Trial.filebase,...
               [diff(sqrt(sum((xyz(:,1,[1,2])-fxyz(:,1,[1,2])).^2,3)));0],...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'ButtWag','bw','b');
nfet.filter('ButFilter',3,5,'low');
nfet.data = [diff(nfet.data);0];

ns = MTADxyz('data',log10(sq(mean(nfet.segs(1:nfet.size(1),50,nan).^2)))','sampleRate',xyz.sampleRate);


eds= linspace(-9,.7,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = ns(ind,1);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = ns(ind,1);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


sts = {'a','w','n','a-w-n'};
hfig = figure;
hfig.Position = [557   274   701   608];
eds= linspace(-9,.7,100);
ads = linspace(-.2,1,100);
for s = 1:4,
    subplot(2,2,s);
    ind = Trial.stc{sts{s}};
    hist2([ns(ind),man(ind)],eds,ads);
    caxis([0,200])
    title(['JPDF: ' ind.label])
    xlabel('Wag_{BL} Power')
    ylabel('PPC_{body} (AU)')
end
saveas(hfig,fullfile('/storage/gravio/figures/req/req20150821',[Trial.filebase '_req20150821_JPDF_WAG_PPC.eps']),'epsc')
saveas(hfig,fullfile('/storage/gravio/figures/req/req20150821',[Trial.filebase '_req20150821_JPDF_WAG_PPC.png']),'png')


sts = {'a','w','n','a-w-n'};
hfig = figure;
hfig.Position = [557   274   701   608];
eds= linspace(-9,.7,100);
ads = linspace(-.8,2,100);
for s = 1:4,
    subplot(2,2,s);
    ind = Trial.stc{sts{s}};
    hist2([ns(ind),fvel(ind)],eds,ads);
    caxis([0,200])
    title(['JPDF: ' ind.label])
    xlabel('Wag_{BL} Power')
    ylabel('Speed_{BL} log10(cm/s)')
end
saveas(hfig,fullfile('/storage/gravio/figures/req/req20150821',[Trial.filebase '_req20150821_JPDF_WAG_SPEED.eps']),'epsc')
saveas(hfig,fullfile('/storage/gravio/figures/req/req20150821',[Trial.filebase '_req20150821_JPDF_WAG_SPEED.png']),'png')
%saveas(hfig,fullfile(Trial.spath,[Trial.filebase '_req20150821_JPDF_WAG_SPEED.eps']),'eps2')


%% REQ #3

Trial = MTATrial('jg05-20120317');


%% Walking features subspace contours
figTitle = 'req20150821_3_JPDF_WAG_PPC_Scontour'

% DEF Mode: 'ss' -> each state has contour 
%           'as' -> primary state has contour, all others states
%                   are merged.
mode = 'ss';


%% Initialize Features
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5,'low');
nfet = MTADfet(Trial.spath,Trial.filebase,...
               [diff(sqrt(sum((xyz(:,1,[1,2])-fxyz(:,1,[1,2])).^2,3)));0],...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'ButtWag','bw','b');
nfet.filter('ButFilter',3,5,'low');
nfet.data = [diff(nfet.data);0];
ns = MTADxyz('data',permute(log10(sq(mean(nfet.segs(1:nfet.size(1), ...
                                                  50,nan).^2))),[2,3,4,1]),'sampleRate',xyz.sampleRate);

man = Trial.load('fet','lsppc');
man.filter('ButFilter',3,1.5,'low');


nbinx = 80;
nbiny = 80;
edx = linspace(-8,-1,nbinx);
edy = linspace(-.2,1,nbiny);

fet = Trial.xyz.copy;
fet.data = [ns(:),man(:)];
% $$$ tag = 'wag_{lower_spine}-vs-speed_{lower_spine}';
% $$$ bhv_JPDF(Trial,ns,fvel,70,70,[.25,1.6],[1.4,2.5],...
% $$$          'WagPow_{spine_lower}',...
% $$$          'Speed_{spine_lower} log10 (cm/s)',{'a-r','r'},...
% $$$          tag);

hfig = figure;
hist2(fet(Trial.stc{'a'},:),edx,edy);
caxis([0,250])

switch mode
  case 'ss'
    sts = {'r','w','s','n','m'};
    stc = 'rwcgm';
  case 'as'
    sts = {'r','a-r'};
    stc = 'rc';
end

hedgs    = {edx};
hedgs(2) = {edy};
edgs    = {edx};
edgs(2) = {edy};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
lbls = {};

hold on,
for i = 1:numel(sts),
    ind = Trial.stc{sts(i)};
    o = hist2(fet(ind,:),hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
    lbls{i} = ind.label;
end
legend(lbls,'Location','NorthWest');

title = 'wag_{lower spine}-vs-PPC_{lower_spine}';
xlabel('WagPow_{spine lower} (AU)');
ylabel('PPC_{body} (AU)');

saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.png']),'png')






% Pitch upperv
%% Rearing feature subspace contours
hostPath = '/gpfs01/sirota/homes/gravio/figures/req/req20150821';
hostPath = '/storage/gravio/figures/req/req20150821';

% TRIAL subset of recording session
Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');
% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);


%% Rear Body pitch vs diff(USpitch)

figTitle = 'req20150821_3_JPDF_USpitch_dUSpitch_Scontour';

% DEF Mode: 'ss' -> each state has contour 
%           'as' -> primary state has contour, all others states
%                   are merged.
mode = 'ss';
clims = [0,150];
contLim = [10,10];


% FET assignment 
fet = Trial.xyz.copy;
fet.data = ([fang(:,3,4,2),[diff(fang(:,3,4,2));0]*fang.sampleRate]);

nbinx = 80;
nbiny = 80;
edx = linspace(-.8,pi/2,nbinx);
edy = linspace(-5,5,nbiny);


hfig = figure;
hist2(fet(Trial.stc{'a'},:),edx,edy);
caxis(clims)

switch mode
  case 'ss'
    sts = {'r','w','s','n','m'};
    stc = 'rwcgm';
  case 'as'
    sts = {'r','a-r'};
    stc = 'rc';
end


hedgs    = {edx};
hedgs(2) = {edy};
edgs    = {edx};
edgs(2) = {edy};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
lbls = {};

hold on,
for i = 1:numel(sts),
    ind = Trial.stc{sts{i}};
    o = hist2(fet(ind,:),hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',contLim,'linewidth',2.5,'Color',stc(i))
    lbls{i} = ind.label;
end
legend(lbls,'Location','NorthWest');
title({'State Contours', ['(Contour Sample Limit: ' num2str(contLim) ')']});
xlabel('pitch_{upper spine} (rad)');
ylabel('d(pitch_{upper spine})/dt (rad/s)');

saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.png']),'png')



%% Rear US pitch vs US log10(abs(diff(ang)))
% REQ 3.2
%
% DEF Mode: 'ss' -> each state has contour 
%           'as' -> primary state has contour, all others states
%                   are merged.

mode = 'ss'; %mode = 'as';
figTitle = 'req20150821_3_JPDF_Bpitch_ladUSpitch_Scontour';
clims = [0,150];
contLim = [10,10];

fet = Trial.xyz.copy;
fet.data = ([fang(:,3,4,2),log10(abs([diff(fang(:,3,4,2));0]*fang.sampleRate))]);


nbinx = 80;
nbiny = 80;
edx = linspace(-.8,pi/2,nbinx);
edy = linspace(-4,.8,nbiny);

hfig = figure;
hist2(fet(Trial.stc{'a'},:),edx,edy);
caxis(clims)

switch mode
  case 'ss';
    sts = {'r','w','s','n','m'};
    stc = 'rwcgm';
  case 'as';  
    sts = {'r','a-r'};
    stc = 'rc';
end


hedgs    = {edx};
hedgs(2) = {edy};
edgs    = {edx};
edgs(2) = {edy};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
lbls = {};

hold on,
for i = 1:numel(sts),
    ind = Trial.stc{sts{i}};
    o = hist2(fet(ind,:),hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',contLim,'linewidth',2.5,'Color',stc(i))
    lbls{i} = ind.label;
end

legend(lbls,'Location','southeast');
title({'State Contours', ['(Contour Sample Limit: ' num2str(contLim) ')']});
xlabel('pitch_{upper spine} (rad)');
ylabel('log10(abs(d(pitch_{upper spine})/dt)) (rad/s)');
h = colorbar;
h.Label.String = ['Sample Count (sr: ' num2str(round(Trial.xyz.sampleRate)) ')'];

saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '_' mode '.png']),'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fet = Trial.xyz.copy;
fet.data = [fang(:,1,2,2),fang(:,3,4,2)];

nbinx = 80;
nbiny = 80;
edx = linspace(.4,pi/2,nbinx);
edy = linspace(-.8,pi/2,nbiny);

hfig = figure;
hist2(fet(Trial.stc{'a'},:),edx,edy);
caxis([0,250])

sts = 'rwsnm';
stc = 'rwcgm';
hedgs    = {edx};
hedgs(2) = {edy};
edgs    = {edx};
edgs(2) = {edy};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
lbls = {};

hold on,
for i = 1:numel(sts),
    ind = Trial.stc{sts(i)};
    o = hist2(fet(ind,:),hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
    lbls{i} = ind.label;
end
legend(lbls,'Location','northwest');

title('State Contours');
xlabel('pitch_{lower spine} (rad)');
ylabel('pitch_{upper spine} (rad)');


saveas(hfig,fullfile(hostPath,[Trial.filebase '_req20150821_3_JPDF_LSpitch_USpitch_Scontour.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_req20150821_3_JPDF_LSpitch_USpitch_Scontour.png']),'png')









hfig = figure; hold on;
plot(fet(Trial.stc{'a'},1),fet(Trial.stc{'a'},2));
hist2(fet(Trial.stc{'a'},:),edx,edy);
cpnts = ClusterPP(hfig);
