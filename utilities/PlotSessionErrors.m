function hfig = PlotSessionErrors(Session)
hfig = figure(8384839);
hfig.Name = 'Transformed distances of ridgid body markers';
hold on
[ep hb hr et] = FindErrorPeriods(Session);
% $$$ subplot(121)
plot([hb.transVec(:,:,2),hr.transVec(:,:,2)])
plot(et)
% $$$ subplot(122)
% $$$ pct = prctile([hb.transVec(:,:,2),hr.transVec(:,:,2)],[1,99]);
% $$$ hist(clip([hb.transVec(:,:,2),hr.transVec(:,:,2)],pct(1),pct(2)),100)

%Lines(ep(:,1),[-300,500],'k');
%Lines(ep(:,2),[-300,500],'r');