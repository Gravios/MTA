function hfig = PlotSessionErrors(Session)
hfig = figure(8384839);
hfig.Name = 'Transformed distances of ridgid body markers';
hold on
[ep hb hr et] = FindErrorPeriods(Session);
plot([hb.transVec(:,:,2),hr.transVec(:,:,2)])
plot(et)
%Lines(ep(:,1),[-300,500],'k');
%Lines(ep(:,2),[-300,500],'r');