function hfig = PlotSessionErrors(Session,varargin)
%function hfig = PlotSessionErrors(Session,hfig,xyz,markers)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('hfig',                             gcf(),                                      ...
                 'xyz',                              'position',                                 ...
                 'markers',   {'head_back','head_left','head_front','head_right'}                ...
);
[hfig,xyz,markers] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

figure(hfig);
if isempty(xyz)&&ischar(xyz),
    xyz = Session.load('xyz',xyz);
end


hfig.Name = 'Transformed distances of ridgid body markers';
hold('on')
[ep hbc et] = FindErrorPeriods(Session,xyz,markers);
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