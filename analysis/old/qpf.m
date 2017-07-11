function pf = qpf(Trial,varargin)
[display,overwrite] = DefaultArgs(varargin,{false,false});

if isa(Trial,'MTATrial'),
pf{1} = MTAPlaceField(Trial,[],{'walk'},overwrite);
pf{2} = MTAPlaceField(Trial,[],{'rear'},overwrite);
pf{3} = MTAPlaceField(Trial,[],{'theta'},overwrite);
pf{4} = MTAPlaceField(Trial,[],{'theta'},overwrite,[],'all','pfcrz');
pf{5} = MTAPlaceField(Trial,[],{'bturn'},overwrite);
pf{6} = MTAPlaceField(Trial,[],{'hturn'},overwrite);
else
pf = Trial;
dispaly=true;
end


if display
    figure,


    while unit~=-1
clf
for i = 1:numel(pf)
        subplotfit(i,numel(pf));
        pf{i}.plot(unit,1);
        title(pf{i}.stateLabel);
end
        unit=figure_controls(gcf,unit);
    end
end

for u = 1:numel(units),
for i = 1:numel(pfsr)
        subplot2(2,numel(pfsr),1,i);
        pfsr{i}.plot(units(u),[],1);
        subplot2(2,numel(pfsc),2,i);
        pfsc{i}.plot(units(u),[],1);
        title(num2str(units(u)));
end
waitforbuttonpress
end


prc = pfsc{4}.plot(94);

prc = 1./log10(abs(prc)).*sign(prc);