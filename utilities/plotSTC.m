function hax = plotSTC(Stc,varargin)
%function hax = plotSTC(Stc,varargin)
%[sampleRate,label_method,states,stateColors] = DefaultArgs(varargin,{Stc.states{1}.sampleRate,'',{},''},true);
[sampleRate,label_method,states,stateColors,staggeredStates] = DefaultArgs(varargin,{Stc.states{1}.sampleRate,'',{},'',true},true);



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

hax = cell([1,nsts]);
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
        hax{i} = patch(xind,repmat([j;j+1;j+1;j],[1,size(xind,2)]),c(i,:));
        hax{i}.EdgeColor = c(i,:);
        hax{i}.FaceAlpha = 1;
        if strcmp(label_method,'text'),
        end
        
    end
end

end

    