function [hfig] = display_session_marker_error_summary(s)

xyz = s.load('xyz');
xper = s.xyz.sync.copy();
xper = xper.data.*s.xyz.sampleRate;
xper = round(xper-xper(1)+1);
xper(end) = xyz.size(1);

% DISPLAY Vicon trial with session errors
hfig = figure();
sp = gobjects([1,0]);
sp(end+1) = subplot(211);
    PlotSessionErrors(s);
    ylim([-100,100]);
sp(end+1) = subplot(212);
    Lines(xper(:,1),[],'g');
    Lines(xper(:,2),[],'r');
    pchIn = gobjects([1,0]);
    pchOut = gobjects([1,0]);
    % GENERATE patches to show mocap recording periods
    for p = 1:size(xper,1),
        pchIn(p) = patch([xper(p,[1,1,2,2])],[-1,1,1,-1],[0.7,0.7,0.7]);
        pchIn(p).FaceAlpha = 0.3;
        if p ~= size(xper,1),
            pchOut(p) = patch([xper(p,[2,2]),xper(p+1,[1,1])],[-1,1,1,-1],[0.1,0.1,0.1]);
            pchOut(p).FaceAlpha = 0.3;
        end
        text(sum(xper(p,:).*[0.8,0.2]),0,{'Vicon Trial',['Ind: ',num2str(p)]});
    end
