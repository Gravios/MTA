
%% REQ #1

Trial = MTATrial('jg05-20120317');

%% Initialize Features
xyz = Trial.load('xyz');
fvel = xyz.vel(1,[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0) = .1;
fvel.data = log10(fvel.data);

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

f = ':';
%eds  = linspace(-.8,2,100);
eds  = linspace(-.2,1,100);
figure, hold on,
fvel = xyz.vel(1,[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0) = .1;
fvel.data = log10(fvel.data);

fet = Trial.xyz.copy;
fet.data = [fvel.data,man.data];


os = [];
ind = Trial.stc{'w'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
hs = bar(eds,histc(os,eds),'histc');
overlay_bar(hs,'b');

os = [];
ind = Trial.stc{'r'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
%hs = bar(eds,histc(os,eds),'histc');
%overlay_bar(hs,'r');

%os = [];
ind = Trial.stc{'n'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
%hs = bar(eds,histc(os,eds),'histc');
%overlay_bar(hs,'g');

%os = [];
ind = Trial.stc{'s'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
%hs = bar(eds,histc(os,eds),'histc');
%overlay_bar(hs,'c');

%os = [];
ind = Trial.stc{'m'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
%hs = bar(eds,histc(os,eds),'histc');
%overlay_bar(hs,'m');


hs = bar(eds,histc(os,eds),'histc');
overlay_bar(hs,'r');




os = [];
ind = Trial.stc{'w'}.data;
for i = ind',os(end+1,:) = median(fet(i',f));end
hs = bar(eds,histc(os,eds),'histc');
overlay_bar(hs,'b');



nfet = MTADfet(Trial.spath,Trial.filebase,...
               [diff(sqrt(sum((xyz(:,1,[1,2])-fxyz(:,1,[1,2])).^2,3)));0],...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'ButtWag','bw','b');
nfet.filter('ButFilter',3,5,'low');
nfet.data = diff(nfet.data);

ns = MTADxyz('data',log10(sq(mean(nfet.segs(1:nfet.size(1),50,nan).^2))),'sampleRate',xyz.sampleRate);


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


