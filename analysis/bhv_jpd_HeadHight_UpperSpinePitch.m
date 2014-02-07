Trial = MTATrial('jg05-20120317');
Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(11)./sum(gausswin(11)));
v = MTADxyz([],[],Filter0(gausswin(31)./sum(gausswin(31)),[zeros(1,9);Trial.vel([],[1,2])]),Trial.xyz.sampleRate);

bins = [50,50]



states  = {'theta','rear&theta','hwalk&theta','lwalk&theta'};


bounds =cat(1,perms([-1.5,1]),[-1.5,-1.5;1,1]);

figure,
for s = 1:numel(states)
subplotfit(s,numel(states))
hist2(cat(1,clip(log10([v(Trial.stc{states{s}},1),v(Trial.stc{states{s}},7)]),-1.5,1),bounds),bins(1),bins(2))
title(['JPD of log10 Speed Lower VS Head Height During ' states{s}]);
end
ForAllSubplots('caxis([0,300])')



legend('rear','high walk','low walk')


xlabel('Sample Count')
ylabel('Height (mm)')

name = 'jpd_lower_spine_head_speed';

saveas(gcf,['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.fig'],'fig')
print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.eps'])
print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/' name '.png'])
