function hax = plotSTC(Stc,varargin)
[sampleRate,label_method] = DefaultArgs(varargin,{Stc.states{1}.sampleRate,''},true);


nsts = numel(Stc.states);
c = jet(nsts);
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

    