Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);

edges = 0:2:100;

states  = {'rear&theta','hwalk&theta','lwalk&theta'};
sc = 'grb';
figure,
bh = [];
for s = 1:numel(states)
    n = histc(Trial.xyz(Trial.stc{states{s},Trial.xyz.sampleRate},1,3),edges);
    hold on
    bh(s) = barh(edges,n,'histc');
    set(bh(s),'FaceColor',sc(s)),set(bh(s),'FaceAlpha',.3)
end


legend('rear','high walk','low walk')

title('Distribution of Lower Spine Height Across Behaviors')
xlabel('Sample Count')
ylabel('Height (mm)')

name = 'lower_spine_height_distrb';

saveas(gcf,['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.fig'],'fig')
print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.eps'])
print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.png'])
