function hax = plotSTC(Stc,varargin)
%function hax = plotSTC(Stc,varargin)
%[sampleRate,label_method,states,stateColors] = DefaultArgs(varargin,{Stc.states{1}.sampleRate,'',{},''},true);
[sampleRate,label_method,states,stateColors] = DefaultArgs(varargin,{Stc.states{1}.sampleRate,'',{},''},true);

Stc = Stc.copy;

if ~isempty(states),
    Stc.states = Stc.states(Stc.gsi(states));
end

nsts = numel(Stc.states);
if isempty(stateColors)||numel(Stc.states)~=numel(stateColors),
    c = jet(nsts);
else
    c = stateColors(:);
end

for i = 1:nsts;
    tper = Stc.states{i}.copy;
    if tper.size(1)>0,
    tper.resample(sampleRate);
    xind = [tper(:,1),tper(:,1),tper(:,2),tper(:,2)]';
    patch(xind,repmat([0;1;1;0],[1,size(xind,2)]),c(i,:));
    if strcmp(label_method,'text'),
    end

    end
end

end

    