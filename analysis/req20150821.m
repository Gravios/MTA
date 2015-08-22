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
legend({'walking trajectory','walking onset','walking offset'})


rtraj = {};
rstart = {};
rstop = {};
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
        pause(.01)
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
