function hax = plotSTC(Stc,varargin)
%function hax = plotSTC(Stc,varargin)
%[sampleRate,label_method,states,stateColors,staggeredStates] =
%DefaultArgs(varargin,{Stc.states{1}.sampleRate,'',{},'',true},true);


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
    states = Stc(states{:});
else
    states = Stc.states;
end

nsts = numel(states);
if isempty(stateColors) || nsts~=numel(stateColors),
    c = jet(nsts);
else
    c = stateColors(:);
end

patchHandle = cell([1,nsts]);
for i = 1:nsts;
    tper = states{i}.copy;
    if staggeredStates,
        j = i;
    else
        j = 0;
    end
    
        
    if tper.size(1)>0,
        tper.resample(sampleRate);
        xind = [tper(:,1),tper(:,1),tper(:,2),tper(:,2)]';
        patchHandle{i} = patch(xind,repmat([j;j+1;j+1;j],[1,size(xind,2)]),c(i,:),'parent',hax);
        patchHandle{i}.EdgeColor = c(i,:);
        patchHandle{i}.FaceAlpha = 1;        
    end
end

hax.YTickMode = 'manual'
hax.YTick = 1.5:numel(states)+0.5;

if strcmp(labelMethod,'text'),
    if staggeredStates,
        set(hax,'YTickLabelMode','manual');
        set(hax,'YTickLabel',cellfun(@(x) x.label,states,'UniformOutput',false));
    end
end

% END MAIN -----------------------------------------------------------------------------------------