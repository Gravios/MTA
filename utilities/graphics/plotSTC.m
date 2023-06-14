function hax = plotSTC(Stc,varargin)
%function hax = plotSTC(Stc,varargin)
%
% DEFARGS :
%    sampleRate,        Stc.states{1}.sampleRate
%    labelMethod,      'text'
%    states,         {{'walk','rear','turn','pause','groom','sit'}}
%    stateColors,      'brgcmy'
%    staggeredStates,   true
%


% DEFARGS ------------------------------------------------------------------------------------------
hax = gca;
defargs = struct(...
    'sampleRate',     Stc.states{1}.sampleRate,...
    'labelMethod',   '',...
    'states',         {{'walk','rear','turn','pause','groom','sit'}},...
    'stateColors',    'brgcmy',...
    'staggeredStates',true ...
);
[sampleRate,labelMethod,states,stateColors,staggeredStates] = ...
    DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

if ~isempty(states),
% ASSUME states is a cell array of state label strings
% COLLECT cell array of MTADepoch objects in the requested order
    states = Stc(states{:});
else
% COLLECT cell array of all MTADepoch objects in MTAStateCollection object
    states = Stc.states;
end

% ASSIGN colors to each state
nsts = numel(states);
if isempty(stateColors)
    c = jet(nsts);
elseif ischar(stateColors)
    c = stateColors(:);
else
    c = stateColors;
end

% CREATE patch objects for each state on the current axis  
patchHandle = cell([1,nsts]);
for s = 1:nsts;
    tper = states{s}.copy;
    j = s.*double(staggeredStates);
        
    if tper.size(1)>0,
        tper.resample(sampleRate);
        xind = [tper(:,1),tper(:,1),tper(:,2),tper(:,2)]';
        patchHandle{s} = patch(xind,repmat([j;j+1;j+1;j],[1,size(xind,2)]),c(s,:),'parent',hax);
        patchHandle{s}.EdgeColor = c(s,:);
        patchHandle{s}.FaceAlpha = 1;        
    end
end


% SET labels along the y-axis associated with each state
hax.YTickMode = 'manual';
hax.YTick = 1.5:numel(states)+0.5;
if strcmp(labelMethod,'text'),
    if staggeredStates,
        set(hax,'YTickLabelMode','manual');
        set(hax,'YTickLabel',cf(@(x) x.label,states));
    end
end


% END MAIN -----------------------------------------------------------------------------------------