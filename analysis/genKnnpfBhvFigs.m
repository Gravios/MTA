sesList = {{'er01-20110719','cof','all'},...
           {'er01-20110721','cof','all'},...
           {'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'walk','hwalk','lwalk'};

numsts = numel(states);
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});

    Pfs = MTAApfs(Trial,[],'walk',1);Pfh = MTAApfs(Trial,[],'hwalk',1);Pfl = MTAApfs(Trial,[],'lwalk',1);
    Pfks = MTAAknnpfs(Trial,[],'walk',1);Pfkh = MTAAknnpfs(Trial,[],'hwalk',1);Pfkl = MTAAknnpfs(Trial,[],'lwalk',1);

    for i = 1:numel(Pfs.data.clu),
        clf
        pfmr = [];
        pkmr = [];

        pfmr(1) = max(Pfs.data.rateMap(:,i,1));
        pfmr(2) = max(Pfh.data.rateMap(:,i,1));
        pfmr(3) = max(Pfl.data.rateMap(:,i,1));

        pkmr(1) = max(Pfks.data.rateMap(:,i,1));
        pkmr(2) = max(Pfkh.data.rateMap(:,i,1));
        pkmr(3) = max(Pfkl.data.rateMap(:,i,1));

        subplot(231);Pfs.plot(i,[],[],[0,max(pfmr)]);title('Walk')
        ylabel('Normal PlaceFields')
        subplot(232);Pfh.plot(i,[],[],[0,max(pfmr)]);title('High Walk')
        subplot(233);Pfl.plot(i,[],[],[0,max(pfmr)]);title('Low Walk')
        subplot(234);Pfks.plot(i,[],[],[0,max(pkmr)]);title('Walk')
        ylabel('Knn PlaceFields')
        subplot(235);Pfkh.plot(i,[],[],[0,max(pkmr)]);title('High Walk')
        subplot(236);Pfkl.plot(i,[],[],[0,max(pkmr)]);title('Low Walk')

            pause(.1),
            name = [strjoin(sesList{ses},'.'), '-walk_pfs_knnpfs-unit-' num2str(i)];
            print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/WalkPlaceFields/' Trial.name '/' name '.eps'])
            print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/WalkPlaceFields/' Trial.name '/' name '.png'])


    end
end