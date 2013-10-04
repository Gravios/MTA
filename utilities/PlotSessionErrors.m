function PlotSessionErrors(Session)
figure
hold on
[ep hb hr et] = FindErrorPeriods(Session);
plot([hb.transVec(:,:,2),hr.transVec(:,:,2)])
plot(et)
%Lines(ep(:,1),[-300,500],'k');
%Lines(ep(:,2),[-300,500],'r');