function disp_ViconTrials(Session)

xyz = Session.xyz.copy;
if xyz.isempty,
    xyz.load(Session);
end

nsync = size(xyz.sync,1);
xyzPeriods = round((xyz.sync.data-xyz.origin)*xyz.sampleRate)+1;

figure,
for i= 1:nsync,
    subplotfit(i,nsync);
    plot(xyz(xyzPeriods(i,1):xyzPeriods(i,2),Session.trackingMarker,1),...
         xyz(xyzPeriods(i,1):xyzPeriods(i,2),Session.trackingMarker,2),...
         '.');
    xlim(Session.maze.boundaries(1,:));
    ylim(Session.maze.boundaries(2,:));
end

         