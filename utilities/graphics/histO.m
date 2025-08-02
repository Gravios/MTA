function [hax,varargout] = histO(varargin)
%function [hax,varargout] = histO(varargin)
% Histograms overlayed
% 
% 
% histO(eds,{{ang(Trial.stc{'a&m'},1,7,3),'m'},...
%            {ang(Trial.stc{'a&s'},1,7,3),'c'},...
%            {ang(Trial.stc{'a&s'},1,7,3),'c'}})
% varargin
%
%   edges
%
%   caX 
%
if isa(varargin{1},'hgsetget'),
    hax = varargin{1};
    varargin(1)=[];
else
    hax = gca;
end

hold(hax,'on');

edges = varargin{1};
caX   = varargin{2};
clear('varargin');

assert(numel(caX)>2,  'MTA:utilities:histO:TooFewInputs')
assert(iscell(caX{1}),'MTA:utilities:histO:IncorrectInputStructure:NotCell')
assert(numel(caX{1})>=2, 'MTA:utilities:histO:IncorrectInputStructure:TooFewInputs')

varargout = cell(1,nargout);
nCaX = numel(caX);

for ind = 1:nCaX,

    % Plot bar graph of each X in caX
    hs = bar(edges,histc(caX{ind}{1},edges),'histc');
    hs.FaceColor = caX{ind}{2};

    % Modify the patch properties of the bar graph
    if numel(caX{ind})>2,
        hs.FaceAlpha = caX{ind}{3};
    else
        hs.FaceAlpha = (1/nCaX+.5)/2;
    end
    
    % Put handles of each bar graph into the variable output
    if nargout>1 && ind<=nargout,
        varargout(ind) = {hs};
    end
    
end

