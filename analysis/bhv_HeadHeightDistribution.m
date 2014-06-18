Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);

edges = 0:10:300;

states  = {'rear&theta','hwalk&theta','lwalk&theta'};
sc = 'grb';
figure,
bh = [];
for s = 1:numel(states)
    n = histc(Trial.xyz(Trial.stc{states{s},Trial.xyz.sampleRate},7,3),edges);
    hold on
    bh(s) = barh(edges,n,'histc');
    set(bh(s),'FaceColor',sc(s)),set(bh(s),'FaceAlpha',.3)
end


legend('rear','high walk','low walk')

title('Distribution of head height across behaviors')
xlabel('Sample Count')
ylabel('Height (mm)')