function hfig = PlotSessionErrors(Session,varargin)
%function hfig = PlotSessionErrors(Session,figureHandle)
%
if ~isempty(varargin),
    hfig = figure(varargin{1})
else    
    hfig = figure(8384839);
end

hfig.Name = 'Transformed distances of ridgid body markers';
hold on
[ep hbc et] = FindErrorPeriods(Session);
% $$$ subplot(121)
plot(et)
for i = 1:numel(hbc)
    plot(hbc{i}.transVec(:,:,2))
end

% $$$ subplot(122)
% $$$ pct = prctile([hb.transVec(:,:,2),hr.transVec(:,:,2)],[1,99]);
% $$$ hist(clip([hb.transVec(:,:,2),hr.transVec(:,:,2)],pct(1),pct(2)),100)

%Lines(ep(:,1),[-300,500],'k');
%Lines(ep(:,2),[-300,500],'r');