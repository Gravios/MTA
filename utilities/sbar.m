function sbar(Trial,fet,varargin)
[states,edgs,hfig,hax] = DefaultArgs(varargin,{{},100,gcf,gca},true);

if iscell(states),
    states = strjoin(states,'-');
end

astates = ['a-' states];

figure(hfig),
hold on
ind = Trial.stc{astates};
bar(hax,edgs,histc(log10(fet(ind,1)),edgs),'histc');
ind = Trial.stc{states};
h = bar(hax,edgs,histc(log10(fet(ind,1)),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;
hold off