function rear_ccg(Session,varargin)
[trialName,ccg_name,display_mode] = DefaultArgs(varargin,{'all','rdpf_nrhp_sur_par20','report'});

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
    Trial = MTATrial(Session,trialName);
end

load([Trial.spath.analysis Trial.filebase '.ccg.' ccg_name '.mat']);

fccg = Bccg.filter(gausswin(7));

numClu = size(Bccg.cluMap,1);

[accg,atbin] = autoccg(Session);
clear('Session');

pf_search = MTAPlaceField([]);
pf_search.mazeName = 'cof';
pf_search.trialName = 'all';
pf_search.trackingMarker = 'head_front';
pf_search.stateLabel = 'head.theta';
pf_search.spk_shuffle = 'n';
pf_search.pos_shuffle = 0;
pf_search.numBSiterations = 1;
pf_search.nbins = 50;
pf_search.smooth = 0.03;

Trial = Trial.load_Pfs();



if ~isempty(display_mode),

    unit = 1;
    figure(2993),
    set(gcf,'CurrentCharacter','l');
    unit = 1;
    while 1,

        %% START - FIG
        if Bccg.partitions>1,
            clf

            subplot2(4,6,1,[1,2])
            imagesc(Bccg.tbin,... 
                    Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,...
                    sq(fccg(:,unit,1,:,1))')
            axis xy
            title(num2str(unit));
            subplot2(4,6,1,[3,4])
            imagesc(Bccg.tbin,...
                    Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,...
                    sq(fccg(:,unit,2,:,1))')
            axis xy
            subplot2(4,6,1,[5,6])
            imagesc(Bccg.tbin,...
                    Bccg.partition_boundaries{3,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{3,1}(unit,:))/2,...
                    sq(fccg(:,unit,3,:,1))')
            axis xy
            subplot2(4,6,2,[1,2])
            mprr = permute(repmat(sq(mean(fccg(10:Bccg.halfBins-10,unit,1,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            imagesc(Bccg.tbin,... 
                    Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,...
                    (sq(fccg(:,unit,1,:,1))./mprr)')
            axis xy
            subplot2(4,6,[3,4],[1,2])
            pf_search.stateLabel = 'rear.theta';
            Pfs = Trial.getPfs(pf_search);
            ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy
            subplot2(4,6,2,[3,4])
            mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,2,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            imagesc(Bccg.tbin,...
                    Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,...
                    (sq(fccg(:,unit,2,:,1))./mprr)')
            axis xy
            subplot2(4,6,2,[5,6])
            mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,3,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            imagesc(Bccg.tbin,...
                    Bccg.partition_boundaries{3,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{3,1}(unit,:))/2,...
                    (sq(fccg(:,unit,3,:,1))./mprr)')
            axis xy
            subplot2(4,6,[3,4],[5,6])
            pf_search.stateLabel = 'head.theta';
            Pfs = Trial.getPfs(pf_search);
            ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy
            subplot2(4,6,[3,4],[3,4])
            bar(atbin,accg(:,unit))
            axis tight
        else

            subplot2(4,6,1,[1,2])
            bar(Bccg.tbin,sq(fccg(:,unit,1,:,1)))
            axis xy
            axis tight
            title(num2str(unit));

            subplot2(4,6,1,[3,4])
            bar(Bccg.tbin,sq(fccg(:,unit,2,:,1)))
            axis tight
            axis xy

            subplot2(4,6,1,[5,6])
            bar(Bccg.tbin,sq(fccg(:,unit,3,:,1)))
            axis tight
            axis xy

            subplot2(4,6,2,[1,2])
            mprr = permute(repmat(sq(mean(fccg(10:Bccg.halfBins-10,unit,1,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            bar(Bccg.tbin,sq(fccg(:,unit,1,:,1))./mprr)
            axis tight
            axis xy

            subplot2(4,6,[3,4],[1,2])
            pf_search.stateLabel = 'rear.theta';
            Pfs = Trial.getPfs(pf_search);
            ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy

            subplot2(4,6,2,[3,4])
            mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,2,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            bar(Bccg.tbin,sq(fccg(:,unit,2,:,1))./mprr)
            axis tight
            axis xy

            subplot2(4,6,2,[5,6])
            mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,3,:,1))),1,size(fccg,1)),[2,1]);
            mprr(mprr==0) = inf;
            bar(Bccg.tbin,sq(fccg(:,unit,3,:,1))./mprr)
            axis tight
            axis xy

            subplot2(4,6,[3,4],[5,6])
            pf_search.stateLabel = 'head.theta';
            Pfs = Trial.getPfs(pf_search);
            ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy
            subplot2(4,6,[3,4],[3,4])
            bar(atbin,accg(:,unit))
            axis tight

        end
        %% END - FIG



        %% START - Controls
        switch display_mode

          case 'report'
            %% ReportFig
            set(gcf,'CurrentCharacter','n');
            reportfig(2993,[Trial.filebase '.ccg.' ccg_name],[],['unit: ' num2str(unit)],[],0);
            unit = unit+1;
            if unit > numClu, return,end

          case 'display'
            %% manual controls
            unit = figure_controls(unit)
            if unit==-1,return,end

        end
        %% END - Controls

    end

end
