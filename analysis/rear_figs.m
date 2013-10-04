function rear_figs(Session,varargin)
[figure_name,figure_type,rear_state_set,trialName,display_mode] = DefaultArgs(varargin,{'','pfs_rearing',{'head1.theta','pre_rear_onset','pre_peri_rear_onset' 'dia_rear_onset','post_peri_rear_onset', 'post_rear_onset','rear.theta';'rear.theta','pre_rear_offset','pre_peri_rear_offset','dia_rear_offset','post_peri_rear_offset','post_rear_offset','head1.theta'},'all','display'});


if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end
if ~isa(Session,'MTATrial'),
    Session = MTATrial(Session,{},trialName);
end


if ~isempty(figure_name),
    figure_name = ['.' figure_name];
end

switch figure_type

  case 'pfs_normal'

    %% Load Place Fields
    Session.Pfs = [];
    Session = Session.load_Pfs;

    pf_search = MTAPlaceField([]);
    pf_search.trackingMarker = 'head_front';
    pf_search.smooth = 0.03;


    %% Display Everything
    figure(102)
    set(gcf,'CurrentCharacter','l');
    unit = 1;
    while 1,
        clf
        for pf = 1:size(rear_state_set,1),
            for s = 1:size(rear_state_set,2)
                pf_search.stateLabel = rear_state_set{pf,s};
                Pfs = Session.getPfs(pf_search);
                %% Rate Map
                subplot2(size(rear_state_set,2),size(rear_state_set,1),s,pf);
                ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit})
                title([Pfs.stateLabel ' ' num2str(unit)])
            end
        end

        numClu = size(Pfs.cluMap,1);


        switch display_mode
          case 'report'
            %% ReportFig
            set(gcf,'CurrentCharacter','n');
            reportfig(102,[Session.filebase '.pfs' figure_name],[],['unit: ' num2str(unit)],[],0);
            unit = unit+1;
            if unit > numClu, return,end
            
          case 'display'
            %% manual controls
            unit = figure_controls(unit);
            if unit==-1,return,end
        end

    end


  case 'pfs_rearing'

    %% Load Place Fields
    Session.Pfs = [];
    Session = Session.load_Pfs;

    pf_search = MTAPlaceField([]);
    pf_search.trackingMarker = 'head_front';
    pf_search.smooth = 0.03;

    event_type = 'onset';
    switch event_type
      case 'onset'
        rear_state_set = {'pre_rear_onset','pre_peri_rear_onset','dia_rear_onset','post_peri_rear_onset','post_rear_onset'};
      case 'offset'
        rear_state_set = {'pre_rear_offset','pre_peri_rear_offset','dia_rear_offset','post_peri_rear_offset','post_rear_offset'};
    end

    figure(102)
    set(gcf,'CurrentCharacter','l');
    unit = 1;
    while 1,
        clf
        for pf = 1:numel(rear_state_set),
            for s = 1:numel(rear_state_set),
                pf_search.stateLabel = rear_state_set{pf};
                Pf1 = Session.getPfs(pf_search);
                pf_search.stateLabel = rear_state_set{s};
                Pf2 = Session.getPfs(pf_search);
                %% Rate Map
                subplot2(numel(rear_state_set),numel(rear_state_set),s,pf);
                ppf(Pf1.xbin,Pf1.ybin,log10((Pf1.rateMap{unit}+1)./(Pf2.rateMap{unit}+1)))
                colorbar
                title([ Pf1.stateLabel(1:end-length(event_type)-6) ...
                        Pf1.stateLabel(end-length(event_type):end) '/' Pf2.stateLabel(1:end-length(event_type)-6) Pf2.stateLabel(end-length(event_type):end) ' ' num2str(unit)],'Interpreter','none')
            end
        end

        numClu = size(Pf2.cluMap,1);
        ForAllSubplots('caxis([-1.3,1.3])');

        switch display_mode
          case 'report'
            %% ReportFig
            set(gcf,'CurrentCharacter','n');
            reportfig(102,[Session.filebase '.pfs_rearing' figure_name],[],['unit: ' num2str(unit)],[],0);
            unit = unit+1;
            if unit > numClu, return,end
            
          case 'display'
            %% manual controls
            unit = figure_controls(unit);
            if unit==-1,return,end
        end

    end



  case 'ccg_rear'

    load([Trial.spath.analysis Trial.filebase '.ccg.' figure_name '.mat']);

    fccg = Bccg.filter(gausswin(7));

    ccg_search = MTAccg([]);
    ccg_search.mazeName = 'cof';
    ccg_search.trialName = 'all';
    ccg_search.name = 'rearing';


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

                % CCG centered on rear onset
                subplot2(4,6,1,[1,2]);imagesc(Bccg.tbin,Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,sq(fccg(:,unit,1,:,1))'),axis xy,title(num2str(unit));
                
                % CCG centered on rear offset
                subplot2(4,6,1,[3,4]);imagesc(Bccg.tbin,Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,sq(fccg(:,unit,2,:,1))'),axis xy
                
                % CCG centered on walking - nrhp
                subplot2(4,6,1,[5,6]);imagesc(Bccg.tbin,Bccg.partition_boundaries{3,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{3,1}(unit,:))/2,sq(fccg(:,unit,3,:,1))'),axis xy

                % CCG centered on rear onset, normalized by average rate between -4 and -2 seconds of ccg
                subplot2(4,6,2,[1,2])
                mprr = permute(repmat(sq(mean(fccg(10:Bccg.halfBins-10,unit,1,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                imagesc(Bccg.tbin,Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,(sq(fccg(:,unit,1,:,1))./mprr)'),axis xy
                
                % CCG centered on rear offset, normalized by average rate between -4 and -2 seconds of ccg
                subplot2(4,6,2,[3,4])
                mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,2,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                imagesc(Bccg.tbin,Bccg.partition_boundaries{2,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{2,1}(unit,:))/2,(sq(fccg(:,unit,2,:,1))./mprr)'),axis xy

                % CCG centered on walking - nrhp, normalized by average rate between -4 and -2 seconds of ccg
                subplot2(4,6,2,[5,6])
                mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,3,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                imagesc(Bccg.tbin,Bccg.partition_boundaries{3,1}(unit,1:Bccg.partitions)+diff(Bccg.partition_boundaries{3,1}(unit,:))/2,(sq(fccg(:,unit,3,:,1))./mprr)'),axis xy
                
                % Place Field: Rearing 
                subplot2(4,6,[3,4],[1,2]); pf_search.stateLabel = 'rear.theta';Pfs = Trial.getPfs(pf_search);
                ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy

                % Place Field: Walk - head
                subplot2(4,6,[3,4],[5,6])
                pf_search.stateLabel = 'head.theta';
                Pfs = Trial.getPfs(pf_search);
                ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy

                % Unit auto ccg
                subplot2(4,6,[3,4],[3,4])
                bar(atbin,accg(:,unit))
                axis tight
            else

                subplot2(4,6,1,[1,2]);bar(Bccg.tbin,sq(fccg(:,unit,1,:,1))),axis xy,axis tight,title(num2str(unit));

                subplot2(4,6,1,[3,4]);bar(Bccg.tbin,sq(fccg(:,unit,2,:,1))),axis tight,axis xy

                subplot2(4,6,1,[5,6]);bar(Bccg.tbin,sq(fccg(:,unit,3,:,1))),axis tight,axis xy

                subplot2(4,6,2,[1,2])
                mprr = permute(repmat(sq(mean(fccg(10:Bccg.halfBins-10,unit,1,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                bar(Bccg.tbin,sq(fccg(:,unit,1,:,1))./mprr),axis tight,axis xy

                subplot2(4,6,2,[3,4])
                mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,2,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                bar(Bccg.tbin,sq(fccg(:,unit,2,:,1))./mprr),axis tight,axis xy

                subplot2(4,6,2,[5,6])
                mprr = permute(repmat(sq(mean(fccg(Bccg.halfBins+10:end-10,unit,3,:,1))),1,size(fccg,1)),[2,1]);
                mprr(mprr==0) = inf;
                bar(Bccg.tbin,sq(fccg(:,unit,3,:,1))./mprr)
                axis tight
                axis xy


                subplot2(4,6,[3,4],[1,2])
                pf_search.stateLabel = 'rear.theta';
                Pfs = Trial.getPfs(pf_search);
                ppf(Pfs.xbin,Pfs.ybin,Pfs.rateMap{unit}),axis xy


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
                reportfig(2993,[Trial.filebase '.ccg.' figure_name],[],['unit: ' num2str(unit)],[],0);
                unit = unit+1;
                if unit > numClu, return,end

              case 'display'
                %% manual controls
                unit = figure_controls(unit);
                if unit==-1;return,end

            end
            %% END - Controls

        end

    end

end
