

sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
statez = {'theta','i_rt','i_wt','i_gt','i_lt'};

numsts = numel(states);
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
genKnnpfBhvFigs.m    for i = 1:numsts,
        [accg,tbin] = autoccg(Trial,[],states{i});

        for u = 1:size(accg,2),
            name = [strjoin(sesList{ses},'.'), '-' statez{i} '-unit-' num2str(u)];
            bar(tbin,accg(:,u));axis tight
            print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/ACCG/' Trial.name '/' name '.eps'])
            print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/ACCG/' Trial.name '/' name '.png'])
        end
    end
end