function MyFirstPaperFigs(Trial,mode,varargin)
% function MyFirstPaperFigs(Trial,mode,varargin)
% Various Figures for the first paper
% 
% for a list of figures and their options enter 'help' into the var mode


switch mode,

  case 'help'
    disp( 'The following figures are available: ')
    disp( '  phase2d   - 2d phase precession ')
    disp( '  phaseTemp - Temporal phase precession ')
    
    
  case 'figure1'
    %% Figure 1 3D tracking -> labeled behavior
    xyz = Trial.load('xyz');
    xyz.filter(gtwin(xyz.sampleRate,.25));
    figure,
    sph = [];
    % A: Image, rat with markers in maze.
    %subplot2(6,6,1,[1,2]); image('rat.png');
    % B: Image, marker skeleton reconstruction
    %subplot2(6,6,1,[3,4]); image('skel.png');
    % C: Image, cameras and maze
    %subplot2(6,6,1,[3,4]); image('skel.png');

    % D: Plot ( Time Series ), COM_head speed ( blue )
    %    Plot ( Time Series ), COM_body speed ( green )
    %    Plot ( Time Series ), COM_head height ( red )

    window_figD = 1300; index_figD = 600;

    coms = cat(2,xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})),...
               xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
    t = linspace(0,size(coms,1)/Trial.xyz.sampleRate,size(coms,1)-1)';
    vcoms = sqrt(sum(diff(coms).^2,3));
    sph(end+1) = subplot2(6,6,2,[1:6]); 
    ax_figD = plotyy(t,clip(coms(1:end-1,2,3),20,350),[t,t],clip(vcoms,0,25));
    set(ax_figD,'xlim',[479,515]);%[1660,1730]);
    set(ax_figD(1),'ylim',[20,280]);
    set(ax_figD(2),'ylim',[-2,10]);



    %% Figure 1 BHV state Place Fields and ccgs
    % rear walk hwalk lwalk

  case 'figure3'
    Trial = MTATrial('jg05-20120309','all');
    Trial.load('nq');
    units = find(Trial.nq.SpkWidthR>0.8&Trial.nq.eDist>18)';
    states = {'rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
    nsts = numel(states);
    for i = 1:nsts,
        pfs{i} =     MTAAknnpfs(Trial,units,states{i},0,'numIter',1000, ...
                                'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
    end


    f = figure(1123);
    for u = pfs{1}.data.clu, 
        for i = 1:nsts,
            subplotfit(i,nsts),pfs{i}.plot(u);title([pfs{i}.parameters.states,':',num2str(u)]),
        end
        pause(.2),
        reportfig(f,[Trial.name,'-pfs_state'],[],['unit: ',num2str(u)]);
    end


    figc = figure(1231);
    for u = pfl.data.clu, 
        for i = 1:nsts,
            subplotfit(i,nsts),pfs{i}.plot(u);title([pfs{i}.parameters.states,':',num2str(u)]),
        end
        pause(.2),
        reportfig(figc,'jg05-201203017-ccg_state',[],['unit: num2str(u)']);
    end

    bccg{i} = gen_bhv_ccg(Trial);
    figure,
    subplot(211),bccg.plot(units(10),1);
    subplot(212),bccg.plot(units(10),2);


    %% End - Figure 1



  case 'bhvJPDF'

    %% Figure 2 - Marker vs Marker Speed Joint Distribution
    MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

    Trial = MTATrial('jg05-20120317','all');
    xyz = Trial.load('xyz');

    states = 'twr';
    states = 'wgl';
    m1 = 'spine_lower';
    m2 = 'head_front';


    nbins=100;mab = 2.1;mib = -2;
    nbins=50;mab = 2.1;mib = 0.5;

    hc = cat(1,mib,mab,mib,mab);
    bc = cat(1,mab,mib,mib,mab);


    dmax = 100;dmin =0.1;
    dxys = xyz.vel([],[1,2]);
    dxys.filter('ButFilter',3,4,'low');

    figure,set(0,'defaulttextinterpreter','none')
    for state=states,subplot(3,1,find(state==states));

        b = clip(log10(dxys(Trial.stc{state},m1)),mib,mab);
        h = clip(log10(dxys(Trial.stc{state},m2)),mib,mab);

        hist2([[b;bc],[h;hc]],nbins,nbins);
        title([m1 ' vs ' m2 ' during ' Trial.stc{state}.label])
        xlabel([m1 ' speed (cm/s)'])
        ylabel([m2 ' speed (cm/s)'])
        caxis([0,250]);
        line([mib,mab],[mib,mab],'color',[1,1,1]);

        xticks = get(gca,'XTickLabel');
        xticks = mat2cell(xticks,ones(1,size(xticks,1)),size(xticks,2));
        xticks = cellfun(@str2num,xticks);
        xticks = mat2cell(10.^xticks,ones(1,size(xticks,1)),1);
        xticks = cellfun(@sprintf,repmat({'%2.2f'},numel(xticks),1),xticks,'uniformoutput',false);
        set(gca,'XTickLabel',xticks);

        yticks = get(gca,'YTickLabel');
        yticks = mat2cell(yticks,ones(1,size(yticks,1)),size(yticks,2));
        yticks = cellfun(@str2num,yticks);
        yticks = mat2cell(10.^yticks,ones(1,size(yticks,1)),1);
        yticks = cellfun(@sprintf,repmat({'%2.2f'},numel(yticks),1),yticks,'uniformoutput',false);
        set(gca,'YTickLabel',yticks);

    end


    %% End - Figure 2


  case 'duno'
    %% Figure 3 - PlaceFields of behaviors for multiple units sorted by behavioral modulation

    sname = 'jg05-20120310.cof.all';
    Trial = MTATrial.validate(sname);
    Trial.load('nq');

    states = {'walk&theta','hwalk&theta','lwalk&theta'};

    % Get the place fields
    pfs = {};
    for s = 1:numel(states);
        pfs{s} = MTAAknnpfs(Trial,[],states{s},0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThresh',70);
    end

    % Get the auto correlograms
    [accg,tbins] = autoccg(MTASession(sname));

    % Select Some units
    spc = [];
    for i=1:numel(pfs),
        spc = cat(2,spc,pfs{i}.spatialCoherence([]));
    end
    mxr = [];
    for i=1:numel(pfs),
        mxr = cat(2,mxr,pfs{i}.maxRate([]));
    end

    ind = (spc(:,2)>.98|spc(:,3)>.98)&Trial.nq.eDist>30&Trial.nq.SpkWidthR>.5;
    figure,plot(log10(mxr(ind,2)),log10(mxr(ind,3)),'.')
    line([-1,2],[-1,2]);

    wunits = find(ind&log10(mxr(:,2))<0);
    wrunits = find(ind&log10(mxr(:,2))>.8&log10(mxr(:,3))>.8);

    [~,wrind] =sort(mxr(wrunits,2)./mxr(wrunits,3));
    wrunits(wrind) = wrunits;

    runits = find(ind&log10(mxr(:,3))<.54);

    sunits = [wunits;wrunits;runits];

    % Plot it all
    nrows = numel(states)+1;
    ncols = numel(sunits);

    figure
    for unit = sunits'
        un = find(unit==sunits);
        mufr = max(mxr(unit,:))
        for s = 1:nrows-1;
            subplot2(nrows,ncols,s,un);
            pfs{s}.plot(unit,[],[],[0,mufr]);
        end
        subplot2(nrows,ncols,s+1,un);
        bar(tbins,accg(:,unit));axis tight;
    end


    %% End Figure - 3

    
  case 'phaseTemp'
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});
    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));

    units = select_units(Trial,18);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.8);

    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;

    states = {'theta','rear&theta','walk&theta','lswalk&theta','hswalk&theta'};
    nsts = numel(states);

    ow = false;
    spk = {};
    for s = 1:nsts,
        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');
        pfs{s} = MTAApfs(Trial,[],states{s},ow);

    end
    spk{1}.create(Trial,xyz.sampleRate,[],[],'deburst');


    hfig = figure(487264);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./100)
    for unit = units,
        clf
        for s = 1:nsts,
            stateOnset = Trial.stc{states{s},xyz.sampleRate}.data(:,1);
            stateOffset = Trial.stc{states{s},xyz.sampleRate}.data(:,2);
            tShiftOnset = round(stateOnset-1*xyz.sampleRate);
            tShiftOffset = round(stateOffset-1*xyz.sampleRate);

            segLen =round(2*xyz.sampleRate);

            myPhaseSegsOnset = GetSegs(tbp_phase(:,phase_chan),tShiftOnset,segLen,0);
            myPhaseSegsOffset = GetSegs(tbp_phase(:,phase_chan),tShiftOffset,segLen,0);


            mySpk = spk{1}(unit); 
            if numel(mySpk) <50,continue,end
            mySpkSegs = false([tbp_phase.size(1),1]);
            mySpkSegs(mySpk) = true;
            
            mySpksegsOnset =  GetSegs(mySpkSegs,tShiftOnset,segLen,false);
            mySpksegsOffset = GetSegs(mySpkSegs,tShiftOffset,segLen,false);
            
            myInd = 1:tbp_phase.size(1);
            myResSegsOnset = GetSegs(myInd,tShiftOnset,segLen,0);
            myResSegsOffset = GetSegs(myInd,tShiftOffset,segLen,0);

            myResSegsOnset = bsxfun(@minus,myResSegsOnset,stateOnset');
            myResSegsOffset = bsxfun(@minus,myResSegsOffset,stateOffset');

            subplot2(3,nsts,1,s)
            pfs{s}.plot(unit,[],1);
            subplot2(3,nsts,2,s)
            try,            hist2([[myResSegsOnset(mySpksegsOnset);myResSegsOnset(mySpksegsOnset)],...
                                   [circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset));circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset))+360],...
                                   [circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset));circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset))+720]],30,25);
% $$$             plot(myResSegsOnset(mySpksegsOnset),circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset)),'.')
% $$$             hold on            
% $$$             plot(myResSegsOnset(mySpksegsOnset),circ_rad2ang(myPhaseSegsOnset(mySpksegsOnset))+360,'.')
% $$$             xlim([-400,400]);
% $$$             ylim([-180,540]);
                subplot2(3,nsts,3,s)
                hist2([[myResSegsOffset(mySpksegsOffset);myResSegsOffset(mySpksegsOffset)],...
                       [circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset));circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset))+360],...
                       [circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset));circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset))+720]],30,25);

% $$$             plot(myResSegsOffset(mySpksegsOffset),circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset)),'.')
% $$$             hold on 
% $$$             plot(myResSegsOffset(mySpksegsOffset),circ_rad2ang(myPhaseSegsOffset(mySpksegsOffset))+360,'.')   
% $$$             xlim([-400,400]);
% $$$             ylim([-180,540]);
            end
        end
        saveas(hfig,['/gpfs01/sirota/home/gravio/figures/bhvPhasePrecession/',[Trial.filebase,'.bpp_temporal-',num2str(unit),'.png']],'png');
        waitforbuttonpress
    end

    
  case 'phaseZ'
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});
    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));
      
    %units =1:100;
    units = select_units(Trial,18);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.8);

    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;

    Trial.stc.states{end+1} = MTAHvel(Trial);
    states = {'theta','vel&theta'};
    %states = {'theta','vel&theta','rear&theta','walk&theta'};    
    nsts = numel(states);

    ow = false;
    spk = {};
    for s = 1:nsts,
        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');
        pfs{s} = MTAApfs(Trial,units,states{s},ow);

    end
    spk{1}.create(Trial,xyz.sampleRate,[],[],'deburst');
    
    %xyz(:,Trial.trackingMarker,3);
    %xyz.filter(gtwin(.2,xyz.sampleRate));
    %ang = Trial.ang.copy;
    %ang.create(Trial,xyz);
    %z = ang(:,'head_back','head_front',2);
    z = xyz(:,Trial.trackingMarker,3);

    aIncr = false;
    hfig = figure(38385);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./50)
    unit = units(1);
    while unit~=-1,
        
        clf
        for s = 1:nsts,
            res = spk{s}(unit);
            
            if numel(res) <50,continue,end
            res(res>numel(z))=[];
            zspk = log10(z(res));
            phzspk = tbp_phase(res,phase_chan);
            
            gind = nniz(zspk)&nniz(phzspk)&isreal(zspk);
            
            subplot2(6,nsts,[1,2],s);
            plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.');
            xlim([-500,500]),ylim([-500,500])
            title(states{s})
            
            subplot2(6,nsts,[3,4],s);
            pfs{s}.plot(unit);
            %hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
            title(num2str(unit))
            
            if sum(gind)>10,
                subplot2(6,nsts,[5,6],s);plot(zspk(gind),circ_rad2ang(phzspk(gind)),'.');
                hold on,          plot(zspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
                %hold on,          plot(zspk(gind),circ_rad2ang(phzspk(gind))+720,'.');
                %xlim([-1.4,1.5]),
                xlim([1.2,2.5]),
                ylim([-180,540])
% $$$                 subplot2(6,nsts,[5,6],s);
% $$$                 try
% $$$                     out = hist2([[zspk(gind);zspk(gind);zspk(gind)],...
% $$$                            [circ_rad2ang(phzspk(gind));...
% $$$                             circ_rad2ang(phzspk(gind))+360;...
% $$$                             circ_rad2ang(phzspk(gind))+720]],linspace(1.4,2.41,30),linspace(-180,900,60));
% $$$                 end
            end
        end
        
        saveas(hfig,['/gpfs01/sirota/home/gravio/figures/bhvPhasePrecession/',[Trial.filebase,'.bpp_1dz-',num2str(unit),'.png']],'png');
        unit = figure_controls(hfig,unit,units,aIncr);
        
        %reportfig(gcf,'er06-20130613-2Dphspredb',0,num2str(unit),[],0);
        
    end


  case 'phase2dDRZ'
    %% Figure 4 - Classical 2-D phase precession
    %Trial = MTATrial(sname,'all');
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});

    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));

    units = select_units(Trial,18);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.8);
    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;


    states = {'theta','rear&theta','walk&theta','lswalk&theta','hswalk&theta'};
    nsts = numel(states);
    ow = true;
    spk = {};
    pfs = {};
    pmr = {};
    pmp = {};
    wpmr = {};
    pfd = {};
    DRZ = {};

    ow = true;

    for s = 1:nsts
        pfs{s} = MTAApfs(Trial,[],states{s},ow,'binDims',[50,50],'SmoothingWeights',[1.5,1.5]);

        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');

        % Get the expected ufr for each xy 
        wpmr{s} = zeros(xyz.size(1),numel(units));
        %wpmr{s} = ones(xyz.size(1),numel(units));
        [~,indx] = min(abs(repmat(pfs{s}.adata.bins{1}',xyz.size(1),1)-repmat(xyz(:,Trial.trackingMarker,1),1,numel(pfs{s}.adata.bins{1}))),[],2);
        [~,indy] = min(abs(repmat(pfs{s}.adata.bins{2}',xyz.size(1),1)-repmat(xyz(:,Trial.trackingMarker,2),1,numel(pfs{s}.adata.bins{2}))),[],2);
        indrm = sub2ind(pfs{s}.adata.binSizes',indx,indy);


        for unit = units,
            rateMap = pfs{s}.plot(unit);
            %rateMap = rot90(rot90(rateMap)');
            wpmr{s}(:,unit==units) = rateMap(indrm);
        end

        % wpmr(wpmr<1)=1;
        % figure,scatter(xyz(1:61:end,7,1),xyz(1:61:end,7,2),wpmr(1:61:end,20))

        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
        pmr{s} = repmat(pmr{s}(:)',xyz.size(1),1);


        %[pr,px] = pfs.maxRate(units)
        %pfs.plot(unit);
        %figure,imagesc(pfs.adata.bins{1},pfs.adata.bins{2},rateMap');
        %hold on,plot(px(unit==units,1),px(unit==units,2),'w*')



        pfds = [];
        for unit = units
            pfhxy = xyz(:,{'head_back','head_front'},:);
            pfhxy = cat(2,pfhxy,permute(repmat([pmp{s}(unit==units,:),0],xyz.size(1),1),[1,3,2]));
            pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
            
            cor = cell(1,3);
            [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
            cor = cell2mat(cor);
            
            por = cell(1,3);
            [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
            por = cell2mat(por);
            
            pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
        end

        pfd{s} = zeros(size(pfds));
        pfd{s}(abs(pfds)<=pi/2)=-1;
        pfd{s}(abs(pfds)>pi/2)=1;


        %DRZ 
        DRZ{s} = pfd{s}.*(1-wpmr{s}./pmr{s});

% $$$         ind = nniz(xyz(:,Trial.trackingMarker,1));
% $$$         figure,
% $$$         plotcc(xyz(ind,Trial.trackingMarker,1),...
% $$$                 xyz(ind,Trial.trackingMarker,2),...
% $$$                 DRZ{1}(ind,11))
% $$$         colorbar
% $$$         xlim([-500,500])
% $$$         ylim([-500,500])

    end

    aIncr = true;
    hfig = figure(38384);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./30)
    unit = units(1);
    while unit~=-1,
        
        clf
        for s = 1:nsts,
            res = spk{s}(unit);
            
            if numel(res) <50,continue,end
            res(res>xyz.size(1))=[];            
            drzspk = DRZ{s}(res,unit==units);
            phzspk = tbp_phase(res,phase_chan);
            
            gind = ~isnan(drzspk)&~isnan(phzspk);
            
            subplot2(6,nsts,[1,2],s);
            plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.');
            xlim([-500,500]),ylim([-500,500])
            title(states{s})
            
            subplot2(6,nsts,[3,4],s);
            pfs{s}.plot(unit);
            hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
            title(num2str(unit))
            
            if sum(gind)>10,
                subplot2(6,nsts,[5,6],s);plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.');
                hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
                hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+720,'.');
                xlim([-1,1]),
                ylim([-180,900])
% $$$                 subplot2(6,nsts,6,s);
% $$$                 hist2([[drzspk(gind);drzspk(gind)],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+720]],30,25);
            end
        end
        
        saveas(hfig,['/gpfs01/sirota/home/gravio/figures/bhvPhasePrecession/',[Trial.filebase,'.bpp_2dDRZ-',num2str(unit),'.png']],'png');
        unit = figure_controls(hfig,unit,units,aIncr);
        
        %reportfig(gcf,'er06-20130613-2Dphspredb',0,num2str(unit),[],0);
        
    end
 

  case 'phase2dRTC'
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});

    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));

    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    units = select_units(Trial,18);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.8);
    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;


    states = {'theta','rear&theta','walk&theta','lswalk&theta','hswalk&theta'};
    nsts = numel(states);
    ow = true;
    spk = {};
    pfs = {};
    pmr = {};
    pmp = {};
    pfd = {};
    pfr = {};
    DRZ = {};

    ow = true;
    for s = 1:nsts,
        pfs{s} = MTAApfs(Trial,[],states{s},ow,'binDims',[50,50],'SmoothingWeights',[1.5,1.5]);

        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');

        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
    end

    units = units(max(cell2mat(pmr),[],2)>5);

    for s = 1:nsts,
        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
    end

    for s = 1:nsts,
        pfrs = [];
        pfds = [];
        for unit = units
            pfhxy = cat(3,xyz(:,{'head_back','head_front'},[1,2]),zeros([xyz.size(1),2]));
            pfhxy = cat(2,pfhxy,permute(repmat([pmp{s}(unit==units,:),0],xyz.size(1),1),[1,3,2]));
            pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
            
            cor = cell(1,3);
            [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
            cor = cell2mat(cor);
            
            por = cell(1,3);
            [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
            por = cell2mat(por);


            [~,~,pfrs(:,unit==units)] = cart2sph(pfhxy(:,3,1)-pfhxy(:,2,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,2,3));


            
            %pfrs(:,unit==units) = por(:,3);
            pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
        end
        
        pfr{s} = pfrs;
        pfd{s} = zeros(size(pfds));
        pfd{s}(abs(pfds)<=pi/2)=-1;
        pfd{s}(abs(pfds)>pi/2)=1;


        %DRZ 
        DRZ{s} = pfd{s}.*pfr{s};

    end


    aIncr = false;
    rtcbins = linspace(-600,600,25);
    hfig = figure(38386);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./30)
    unit = units(1);
    while unit~=-1,
        
        clf
        for s = 1:nsts,
            res = spk{s}(unit);
            
            if numel(res) <50,continue,end
            res(res>xyz.size(1))=[];            
            drzspk = DRZ{s}(res,unit==units);
            phzspk = tbp_phase(res,phase_chan);
            
            gind = ~isnan(drzspk)&~isnan(phzspk);
            
            subplot2(6,nsts,[1,2],s);
            plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.');
            hold on 
            quiver(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),cos(ang(res,'head_back','head_front',1)),sin(ang(res,'head_back','head_front',1)));
            xlim([-500,500]),ylim([-500,500])
            title(states{s})
            
            subplot2(6,nsts,[3,4],s);
            pfs{s}.plot(unit);
            hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
            title(num2str(unit))
            

                subplot2(6,nsts,[5,6],s);
% $$$                 plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.');
% $$$                 hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
% $$$                 hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+720,'.');
% $$$                 %xlim([-1,1]),
% $$$                 ylim([-180,900])
% $$$                 subplot2(6,nsts,6,s);
% $$$                 hist2([[drzspk(gind);drzspk(gind)],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+720]],30,25);
                
                [~,rbin] = histc(drzspk,rtcbins);
                gind = gind&rbin~=0;
            if sum(gind)>10,
                phm = circ_rad2ang(accumarray(rbin(gind),phzspk(gind),size(rtcbins'),@circ_mean));
                phs = circ_rad2ang(accumarray(rbin(gind),phzspk(gind),size(rtcbins'),@circ_std));
                errorbar(rtcbins,phm,phs,phs,'.')
                hold on
                errorbar(rtcbins,phm+360,phs,phs,'.')
                
            end
        end
        
        saveas(hfig,['/gpfs01/sirota/home/gravio/figures/bhvPhasePrecession/',[Trial.filebase,'.bpp_2dRTCns-',num2str(unit),'.png']],'png');
        unit = figure_controls(hfig,unit,units,aIncr);
        
        %reportfig(gcf,'er06-20130613-2Dphspredb',0,num2str(unit),[],0);
        
    end
 

  case 'phase3dRTC'
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});

    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));

    units = select_units(Trial,18);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.8);
    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;


    states = {'theta','rear&theta','walk&theta','lswalk&theta','hswalk&theta'};
    nsts = numel(states);
    ow = true;
    spk = {};
    pfs = {};
    pmr = {};
    pmp = {};
    pfd = {};
    pfr = {};
    pfp = {};
    DRZ = {};

    ow = true;
    for s = 1:nsts,
        pfs{s} = MTAApfs(Trial,[],states{s},ow,'binDims',[50,50,50],'SmoothingWeights',[1.5,1.5,1.5],'type','xyz');

        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');

        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
    end

    units = units(max(cell2mat(pmr),[],2)>5);

    for s = 1:nsts,
        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
    end

    for s = 1:nsts,
        pfrs = [];
        pfps = [];
        pfds = [];
        for unit = units
            pfhxy = xyz(:,{'head_back','head_front'},:);
            pfhxy = cat(2,pfhxy,permute(repmat([pmp{s}(unit==units,:)],xyz.size(1),1),[1,3,2]));
            pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
            
            cor = cell(1,3);
            [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
            cor = cell2mat(cor);
            
            por = cell(1,3);
            [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
            por = cell2mat(por);

            pfrs(:,unit==units) = por(:,3);
            pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
            pfps(:,unit==units) = circ_dist(cor(:,2),por(:,2));        
        end
        
        pfr{s} = pfrs;
        tangents = cat(3,sin(pfds).*cos(pfps), sin(pfds).*sin(pfps),cos(pfds));
        mxyz = repmat(pfhxy(:,2,:)-pfhxy(:,1,:),[1,size(tangents,2),1]);
        pfds = acos(dot(tangents,mxyz,3)./(sqrt(sum(tangents.^2,3)).*sqrt(sum(mxyz.^2,3))));

        pfd{s} = zeros(size(pfds));
        pfd{s}(abs(pfds)<=pi/2)=1;
        pfd{s}(abs(pfds)>pi/2)=-1;


        %DRZ 
        DRZ{s} = pfd{s}.*pfr{s};

    end

    aIncr = true;
    hfig = figure(38387);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./30)
    unit = units(1);
    while unit~=-1,
        
        clf
        for s = 1:nsts,
            res = spk{s}(unit);
            
            if numel(res) <50,continue,end
            res(res>xyz.size(1))=[];            
            drzspk = DRZ{s}(res,unit==units);
            phzspk = tbp_phase(res,phase_chan);
            
            gind = ~isnan(drzspk)&~isnan(phzspk);
            
            subplot2(6,nsts,[1,2],s);
            plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.');
            xlim([-500,500]),ylim([-500,500])
            title(states{s})
            
            subplot2(6,nsts,[3,4],s);
            pfs{s}.plot(unit,'xy');
            hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
            title(num2str(unit))
            
            if sum(gind)>10,
                subplot2(6,nsts,[5,6],s);plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.');
                hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
                hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+720,'.');
                %xlim([-1,1]),
                ylim([-180,900])
% $$$                 subplot2(6,nsts,6,s);
% $$$                 hist2([[drzspk(gind);drzspk(gind)],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
% $$$                        [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+720]],30,25);
            end
        end
        
        saveas(hfig,['/gpfs01/sirota/home/gravio/figures/bhvPhasePrecession/',[Trial.filebase,'.bpp_3dRTC-',num2str(unit),'.png']],'png');
        unit = figure_controls(hfig,unit,units,aIncr);
        
        %reportfig(gcf,'er06-20130613-2Dphspredb',0,num2str(unit),[],0);
        
    end
 

  case 'phase3dRTCVTC'
    [chans,phase_chan] = DefaultArgs(varargin,{4,1});

         
    %vars for er06-20130612
    Trial = MTATrial('er06-20130612');
    units = [33,87,115,121,151];
    units = [115];
    chans = 4; phase_chan = 1;
    
    Trial = MTATrial('er06-20130614','all-cof');
    chans = 4; phase_chan = 1;
    
     %vars for jg05-20120310
    Trial = MTATrial('jg05-20120310');
    units = [10,13,20,25,42,69];
    chans = [72]; phase_chan = 1;

    %vars for jg05-20120310
    Trial = MTATrial('jg05-20120309');
    %units = [10,13,20,25,42,69];
    chans = [72]; phase_chan = 1;


    units = select_units(Trial,20);
    Trial.load('nq');
    units = units(Trial.nq.SNR(units)>.6);

    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.25,xyz.sampleRate));
    
    lfp = Trial.lfp.copy;
    lfp.load(Trial,chans);
    lfp.resample(xyz);
    tbp_phase = lfp.phase;

    ow = true;    
    %states = {'theta','rear&theta','walk&theta','lswalk&theta','hswalk&theta'};
    states = {'walk'};
    states = {'theta'};
    nsts = numel(states);

    spk = {};    
    pfs = {};    pmr = {};    pmp = {};
    DRZ = {};    VTC = {};

    for s = 1:nsts,
        pfs{s} = MTAApfs(Trial,[],'theta',ow,'binDims',[50,50,50],'SmoothingWeights',[1.5,1.5,1.5],'type','xyz');
        %pfs{s} = MTAApfs(Trial,[],states{s},ow,'binDims',[50,50],'SmoothingWeights',[1.5,1.5],'type','xy');
        [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
        pmp{s}(pmp{s}(:,3)<50,3) = 75;        

        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');
    end

%    units = units(max(cell2mat(pmr),[],2)>5);
%     for s = 1:nsts,
%         [pmr{s},pmp{s}] = pfs{s}.maxRate(units);
%     end

    %get DRZ and VTC
    drzdims = [1,2,3];
    txyz = xyz.copy;
    txyz.data = xyz(:,{'head_back','head_front'},drzdims);
    sxyz = xyz.copy;
    sxyz.data = cat(2,xyz(:,'head_back',drzdims),circshift(xyz(:,'head_back',drzdims),-round(xyz.sampleRate/2)));
    for s = 1:nsts,
        %DRZ{s} = pfDRZ(txyz,pmp{s});
        %VTC{s} = pfDRZ(sxyz,pmp{s});
        DRZ{s} = pfDRD(txyz,pmp{s});
        VTC{s} = pfDRD(sxyz,pmp{s});
    end
    
    vel = xyz.vel('head_front',[1,2]);
    vel.data = log10(vel.data);

    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    ael = circ_dist(ang(:,'head_back','head_front',1),circshift(ang(:,'head_back','head_front',1),10));
    ael = ang(:,'head_back','head_front',2);

    [accg,tbins] = autoccg(Trial,units);


    %aIncr = true;
    phase_chan = 1;
    aIncr = false;
    hfig = figure(38338);
    set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]./30)
    unit = units(1);
    pname = '3dVTDZ';
    %unit = 71;
    while unit~=-1,
        
        clf
        for s = 1:nsts,
            res = spk{s}(unit);
            
            if numel(res) <150,continue,end
            res(res>xyz.size(1))=[];            
            drzspk = DRZ{s}(res,unit==units);
            vtcspk = VTC{s}(res,unit==units);
            zspk = xyz(res,'head_front',3);
            vspk = vel(res);
            aspk = ael(res);
            phzspk = tbp_phase(res,phase_chan);
            
            gind = ~isnan(drzspk)&~isnan(vtcspk)&~isnan(phzspk);
            
            subplot2(3,3,1,1);
            plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.');
            xlim([-500,500]),ylim([-500,500])
            title(states{s})
            
            subplot2(3,3,2,1);
            pfs{s}.plot(unit,'xy');
            hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
            title(num2str(unit))
            
            subplot2(3,3,3,1);
            bar(tbins,accg(:,unit)),
            axis tight
            title(['accg: ' num2str(unit)]);
            if sum(gind)>10,                
                %chsv = jet;
                chsv = hsv;
                cim = linspace(-pi,pi,64);
                [xx,xi] = NearestNeighbour(cim,phzspk(gind));
                
                %% 2D VTC of head direction
                subplot2(3,3,1,2);
                %x = log10(abs(drzspk(gind))-2).*sign(drzspk(gind));     xl=[-.7,.7];
                x = drzspk(gind);                                     xl=[-500,500];
                scatter(x,log10(zspk(gind)),7,chsv(xi,:),'filled');xlim(xl)
                ylim([1.5,2.5])
                title('phase(drzH,log10(Z))')
                subplot2(3,3,2,2);

                drzbins = -500:20:500;
                zbins = 1.4:.04:2.5;
                [~,rbin] = histc(drzspk,drzbins);
                [~,zbin] = histc(log10(zspk),zbins);
                tind = (gind&rbin~=0)&(gind&zbin~=0);
                A = accumarray([rbin(tind),zbin(tind)],phzspk(tind),[numel(drzbins),numel(zbins)],@circ_mean,nan);
                imagescnan({drzbins,zbins,A'},[],1,1,[0,0,0]);axis xy
                title('mean phase(drzH,log10(Z))')
                
                %scatter(drzspk(gind),log10(zspk(gind)),7,chsv(xi,:),'filled');
                %xlim([-500,500]),ylim([1.4,2.5])
                subplot2(3,3,3,2);
                hist2([[drzspk(gind);drzspk(gind)],...
                      [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360]],40,30);
                title('JPDF drzH vs phase')
%                 plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.'); hold on
%                    plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
%                    ylim([-180,540])

                %% 2D VTC of Trajectory direction
                subplot2(3,3,1,3);
                %x = log10(abs(vtcspk(gind))-2).*sign(vtcspk(gind));  xl=[-.7,.7];
                x = vtcspk(gind);                                  xl=[-500,500];
                scatter(x,log10(zspk(gind)),7,chsv(xi,:),'filled');xlim(xl);
                ylim([1.5,2.5])
                title('phase(drzT,log10(Z))')
                subplot2(3,3,2,3);
                %scatter(vtcspk(gind),log10(zspk(gind)),7,chsv(xi,:),'filled');
                
                vtcbins = -500:20:500;
                zbins = 1.4:.04:2.5;
                [~,rbin] = histc(vtcspk,vtcbins);
                [~,zbin] = histc(log10(zspk),zbins);
                tind = (gind&rbin~=0)&(gind&zbin~=0);
                A = accumarray([rbin(tind),zbin(tind)],phzspk(tind),[numel(vtcbins),numel(zbins)],@circ_mean,nan);
                imagescnan({vtcbins,zbins,A'},[],1,1,[0,0,0]);axis xy
                title('mean phase(drzT,log10(Z))')
                %xlim([-500,500]),ylim([1.4,2.5])
                subplot2(3,3,3,3);
                hist2([[vtcspk(gind);vtcspk(gind)],...
                      [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360]],40,30);
                title('JPDF drzT vs phase')
%                 plot(vtcspk(gind),circ_rad2ang(phzspk(gind)),'.'); hold on
%                    plot(vtcspk(gind),circ_rad2ang(phzspk(gind))+360,'.');
%                    ylim([-180,540])

            end
        end
        

        figname = ['/gpfs01/sirota/home/gravio/figures/PIM20150112/PhasePrecession/',[Trial.filebase,'.bpp_',pname,'-',num2str(unit),'.png']];
        unit = figure_controls(hfig,unit,units,aIncr,figname);
        
    end 

    

    ufr = Trial.ufr.copy;
    ufr.create(Trial,xyz,'theta',units,.25);

    figure,
    sp(1) = subplot(211);
    plotyy(1:xyz.size-1,diff(Filter0(gtwin(5,xyz.sampleRate),spor(:,3))),1:xyz.size,por(:,3)),ylim([-5,5]),Lines([],0,'k');Lines([],-2.5,'k');
    hold on,Lines(Trial.stc{'w'}(:),[],'k');
    hold on,Lines(Trial.stc{'r'}(:),[],'r');
    Lines(res(phzspk<0&drzspk>100),[],'m');
    sp(2) = subplot(212);
    plot([res;res],[phzspk;phzspk+2*pi]*180/pi,'.');
    linkaxes(sp,'x');
    

    

    figure,sp = [];
    sp(1) = subplot(211);
    plot(find(DRZ{1}(:,1)<0),DRZ{1}(DRZ{1}(:,1)<0,1),'.r')
    hold on
    plot(find(DRZ{1}(:,1)>0),-DRZ{1}(DRZ{1}(:,1)>0,1),'.b')
    sp(2) = subplot(212),
    plot(find(DRZ{1}(:,1)<0),DRZ{1}(DRZ{1}(:,1)<0,1),'.r')
    hold on
    plot(find(DRZ{1}(:,1)>0),-DRZ{1}(DRZ{1}(:,1)>0,1),'.b')
    linkaxes(sp,'xy');
    
    %% break it down to the trajectories
    unit = 115;
    mar = 'head_front';

    res = spk{1}(unit);    
    sper = Trial.stc{'w',xyz.sampleRate}.copy;
    sper.data = [res-1,res+1];
    sper = sper+[-.5,.5];
    sper.cast('TimeSeries');
      
   [pmr{1},pmp{1}] = pfs{s}.maxRate(units,true);
    pfc = pmp{1}(units==unit,:);
    pfcd = sqrt(sum(bsxfun(@minus,sq(xyz(:,mar,:)),pfc).^2,2));
    
    ufr = Trial.ufr.copy;
    ufr.create(Trial,xyz,'theta',115,.5);

    
    figure,
    
    sp(1) = subplot(1,2,1);
    ind = sign(DRZ{1}(:,units==unit))==-1&pfcd<300&sper&xyz(:,7,3)<180;
    c = jet(100);
    %scatter3(xyz(ind,mar,1),xyz(ind,mar,2),xyz(ind,mar,3),7,c(ufr(ind)+1,:))
    scatter(xyz(ind,mar,1),xyz(ind,mar,2),7,c(ufr(ind)+1,:))
    hold on,
    scatter(pfc(1),pfc(2),pfc(3));
    
    sp(2) = subplot(1,2,2);
    ind = sign(DRZ{1}(:,units==unit))==1&pfcd<300&sper&xyz(:,7,3)<180;
    %scatter3(xyz(ind,mar,1),xyz(ind,mar,2),xyz(ind,mar,3),7,c(ufr(ind)+1,:))
    scatter(xyz(ind,mar,1),xyz(ind,mar,2),7,c(ufr(ind)+1,:))
    hold on,
    scatter(pfc(1),pfc(2),pfc(3));
    linkprop(sp,{'cameraposition','cameraupvector'});


    
    %% End Figure 4


  case 'markerDistJPDF'
    xyz = Trial.load('xyz');
    xyz.filter(gtwin(.1,xyz.sampleRate));
    ang = Trial.ang.copy;
    ang.create(Trial,xyz);

    figure,plot(diff(ang(:,'spine_lower','pelvis_root',3)-ang(:,'pelvis_root','spine_middle',3)))
    hold on,plot(diff(ang(:,'pelvis_root','spine_middle',3)-ang(:,'spine_lower','spine_upper',3)),'r')
    figure,plot(diff(ang(:,'spine_lower','pelvis_root',3)-ang(:,'pelvis_root','spine_middle',3))-diff(ang(:,'pelvis_root','spine_middle',3)-ang(:,'spine_lower','spine_upper',3)))

    hold on,Lines(Trial.stc{'w'}(:),[],'k');
    
    fet = MTADxyz('data',[diff([0;Filter0(gausswin(5)./sum(gausswin(5)),diff(ang(:,'pelvis_root','spine_middle',3)-ang(:,'spine_lower','spine_middle',3)))]);0],'sampleRate',ang.sampleRate);
    fet = MTADxyz('data',[diff([0;Filter0(gausswin(5)./sum(gausswin(5)),diff(ang(:,'pelvis_root','spine_middle',2)-ang(:,'spine_lower','spine_middle',2)))]);0],'sampleRate',ang.sampleRate);
    
    [ys,fs,ts] = fet_spec(Trial,fet,[],'wcsd','overwrite',true);
    figure,imagesc(ts,fs,nunity(log10(ys.data))'),axis xy

    [U,S,V] = svd(cov(log10(ys(nniz(ys),:))));

    nfet = zeros([ys.size(1),1]);
    nfet(nniz(ys)) = log10(ys(nniz(ys),:))*V(:,1);
    pfet1 = nfet;
    pfet2 = nfet;
    figure,plot(nfet)
    hold on,Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k');

    figure,hold on
    plot([diff([0;Filter0(gausswin(5)./sum(gausswin(5)),diff(ang(:,'pelvis_root','spine_middle',3)-ang(:,'spine_lower','spine_middle',3)))]);0],'b')
    plot([diff([0;Filter0(gausswin(5)./sum(gausswin(5)),diff(ang(:,'pelvis_root','spine_middle',2)-ang(:,'spine_lower','spine_middle',2)))]);0],'r')



    figure,
    mp = {{'spine_lower','pelvis_root'},{'pelvis_root','spine_upper'};...
          {'spine_lower','pelvis_root'},{'spine_middle','spine_upper'};...
          {'spine_middle','pelvis_root'},{'spine_lower','spine_upper'};...
          {'spine_middle','spine_upper'},{'pelvis_root','spine_upper'}};
    ind = {'gper','walk','rear'};
    for i = 1:size(mp,1),
        for j = 1:numel(ind),
        subplot2(numel(ind),size(mp,1),j,i);
        if j == 1,
            U = cell(3,1);
            Am = [];
            As = [];            
        else
            U = cell(1,1); 
        end
        [U{:}] = nunity([ang(Trial.stc{ind{j}},...
                          mp{i,1}{1},...
                          mp{i,1}{2},3),...
                      ang(Trial.stc{ind{j}},...
                          mp{i,2}{1},...
                          mp{i,2}{2},3)],...
                     @inf,Am,As);
        if j==1,
            Am = U{2};
            As = U{3}; 
        end
        hist2(U{1},-2:.05:2,-2:.05:2);
        end
    end
    

    
    %moving PCA maximization of breathing feature
    [U,mu,vars] = pca(cov(sq(xyz(912800:914000,3,:))));
    txyz = multiprod(U,sq(xyz(:,1,:)),[1,2],2);

  case 'StateUFR'
    %% Figure 5 - State Wise Unit Firing Rates
    sname = 'jg05-20120309';
    Trial = MTATrial(sname,'all');
    Trial.load('nq');

    states = {'theta','vel&theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
    units = find(Trial.nq.eDist>30&Trial.nq.SpkWidthR>.5)';
    sscount = nan(numel(units),numel(states));
    ssdur   = nan(numel(units),numel(states));

    for s = 1:numel(states),

        tsts = Trial.stc{states{s}};
        Trial.spk.create(Trial,Trial.xyz.sampleRate,states{s},units);
        ssdur(:,s) = sum(diff(tsts.data,1,2));
        for u = units,

            try, sscount(u==units,s) = numel(Trial.spk(u));end

        end
    end
    srates = sscount./(ssdur./Trial.xyz.sampleRate);

    nrates = srates./repmat(max(srates,[],2),1,numel(states));

    [~,rind] = sort(nrates(:,3));

    figure,imagescnan(nrates(rind,:)',[],[],1);
    set(gca,'YTickLabel',states);

    %% End Figure 5



  case 'lfp_psd_jpdf'
    %% Figure 6 - 

    Trial = 'jg05-20120317.cof.all';
    %sname = 'jg05-20120310';
    %sname = 'jg05-20120309';
    display = 0;
    %chans = [1:2:8];
    chans = [71:2:96];
    marker = 'spine_lower'

    Trial = MTATrial.validate(Trial); 

    xyz = Trial.load('xyz');
    ang = create(MTADang,Trial,xyz);
    lfp = Trial.load(lfp,chans);


    wlfp = WhitenSignal(lfp.data,[],1);


    pobj = parpool(8);

    tl=[];
    fl=[];
    yl=[];
    parfor i = 1:lfp.size(2),
        [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),2^12,lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
    end

    th=[];
    fh=[];
    yh=[];
    parfor i = 1:lfp.size(2),
        [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),2^9,lfp.sampleRate,2^8,2^8*0.875,[],[],[],[40,120]);
    end
    yld = MTADlfp('data',yl,'sampleRate',1/diff(tl(1:2,1)));
    yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2,1)));

    bang = ButFilter(ang(:,4,5,3),3,[1,20]./(ang.sampleRate./2),'bandpass');
    bhh = ButFilter(xyz(:,7,3),3,[1,20]./(xyz.sampleRate./2),'bandpass');
    bhx = ButFilter(xyz(:,7,1),3,[1,20]./(xyz.sampleRate./2),'bandpass');
    bhy = ButFilter(xyz(:,7,2),3,[1,20]./(xyz.sampleRate./2),'bandpass');
    if display, figure,plot([bang,bhh,bhx,bhy]+3.*repmat(1:4,size(bang,1),1)), end

    wang = WhitenSignal([bang,bhh,bhx,bhy],[],1);
    [ya,fa,ta,phia,fsta] = mtchglong(wang,2^9,ang.sampleRate,2^8,2^8*0.875,[],[],[],[2,16]);

    %figure,plot(mean(ya(:,fa>9&fa<12),2)./mean(ya(:,fa<7),2))
    %figure,plot([mean(ya(:,fa>9&fa<12,1,1),2),mean(ya(:,fa>9&fa<12,2,2),2)])
    %figure,plot(mean(ya(:,fa>9&fa<12,1,1),2),mean([mean(ya(:,fa>9&fa<12,4,4),2),mean(ya(:,fa>9&fa<12,3,3),2)],2),'.')
    %figure,plot(log10(mean(ya(:,fa>9&fa<12,1,1),2)),log10(mean([mean(ya(:,fa>9&fa<12,4,4),2),mean(ya(:,fa>9&fa<12,3,3),2)],2)),'.')
    spowa = log10(mean(ya(:,fa>6&fa<12,1,1),2));
    spowa = log10(mean(ya(:,fa>6&fa<12,2,2),2));



    % VELOCITY 
    xyz = xyz.copy;
    xyz.filter(gausswin(31)./sum(gausswin(31)));
    v = MTADxyz([],[],sqrt(sum(diff(xyz(:,marker,[1,2])).^2,3)).*xyz.sampleRate./10,xyz.sampleRate);



    %figure,hist2([clip(log10(v.data),-2,3),clip(spow,-6,1)],100,100),caxis([0,40])
    %figure,hist2([clip(log10(v.data),-2,3),clip(mean([spowa,spowh],2),-6,1)],100,100),caxis([0,40])

    yad = MTADlfp([],[],ya,1/diff(ta(1:2)));

    % RESAMPLE variables
    xyz.resample(yad);
    v.resample(yad);
    yld.resample(yad);
    %hdv.resample(yad);


    sbins = 25;
    sedges = [-5,-2];
    %sedges = [50,160];

    vbins = 25;
    vedges = [0,2];

    wper = Trial.stc{'w'}.copy;
    wper.cast('TimeSeries');
    wper.resample(yad);

    %spow = xyz(:,7,3);
    spow = spowa;
    spow = clip(spow,sedges(1),sedges(2));
    sind = spow>sedges(1)&spow<sedges(2);
    vlog = clip(log10(v.data),vedges(1),vedges(2));
    vind = vlog>vedges(1)&vlog<vedges(2)&~isnan(vlog);
    aind = sind&vind&wper.data;

    [spow_count,shind] = histc(spow(aind),sedges(1):abs(diff(sedges))/sbins:sedges(2));
    [v_count,vhind] = histc(vlog(aind),vedges(1):abs(diff(vedges))/vbins:vedges(2));

    %tpow = log10(mean(yld(aind,fl>=6&fl<=12,3),2));
    tpow = log10(mean(yld(aind,fl>=4&fl<=16,1,1),2));
    %tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
    tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;

    %A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
    %figure,
    %subplot(133),imagescnan({vedges,sedges,A'},[-5,-3],[],1,[0,0,0]),axis xy,
    %clf
    %imagescnan({vedges,sedges,A'},[1.2,1.9],[],1,[0,0,0]),axis xy,
    %subplot(122),imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,

    % CA1 LM 81;
    % DG  G  85;
    % CA3 ?  95;
    chan = find(chans == 71);
    numIter = 10000;
    %tpow = log10(mean(yld(aind,fl>6&fl<12,chan),2));
    %tpow = log10(mean(yld(aind,fl<=4,chan),2));
    tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
    %tpow = log10(mean(yld(aind,fh>50&fh<80,chan),2));

    B=nan(vbins,sbins,numIter);
    A=nan(vbins,sbins,numIter);
    %S=nan(vbins,sbins,numIter);
    tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
    B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
    A(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
    %S(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
    for i = 2:numIter,
        tpow = tpow(randperm(numel(tpow)));
        tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
        A(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
        %S(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
    end

    AS = sort(A,3);
    P = 1./sum(repmat(A(:,:,1),[1,1,numIter])>A,3);
    P(isinf(P)) = nan;

    SIG = P<=0.0002;
    ASIG = A; 
    ASIG(~SIG)=nan;
    ASIG(B<10)=nan;
    Aclims = [prctile(ASIG(~isnan(ASIG)),5),prctile(ASIG(~isnan(ASIG)),95)];

    figure

    subplot(131),imagescnan({vedges,sppedges,B'./yad.sampleRate},[],[],1,[0,0,0]),axis xy,
    title('Occupancy in seconds')
    ylabel('log10(P(10Hz)) of Spine to Head Distance')
    xlabel('Head Speed (cm/s)')
    ticks_lin2log(gca,'x')

    subplot(132),imagescnan({vedges,sedges,   A(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
    title('Mean Power 1-4Hz given 10Hz Osc. Power dist(SU,HB) VS Vel(HF)')
    ylabel('log10(P(10Hz)) of Spine to Head Distance')
    xlabel('Head Speed (cm/s)')
    ticks_lin2log(gca,'x')

    subplot(133),imagescnan({vedges,sedges,ASIG(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
    title('P<0.0002 and Occupancy > 1.33 seconds')
    ylabel('log10(P(10Hz)) of Spine to Head Distance')
    xlabel('Head Speed (cm/s)')
    ticks_lin2log(gca,'x')





    flim=[1,4;6,12;20,27;30,40];
    %flim=[40,60;60,80;80,100;100,120];
    %mychans = [71,73,81,85,95];
    mychans = 1:2:8;
    B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
    figure
    for c = 1:numel(mychans),
        for i = flim',
            subplot2(numel(mychans),size(flim,1),c,find(i(1)==flim(:,1)));
            tpow = log10(mean(yld(aind,fl>i(1)&fl<i(2),find(chans==mychans(c))),2));
            %tpow = log10(mean(yld(aind,fh>i(1)&fh<i(2),find(chans==mychans(c))),2));
            tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
            AFB = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
            AFB(B<5)=nan;
            AFBclims = [prctile(AFB(~isnan(AFB)),5),prctile(AFB(~isnan(AFB)),95)];
            imagescnan({vedges,sedges,AFB'},AFBclims,[],1,[0,0,0]);
            axis xy,
            ticks_lin2log(gca,'x')
            title(['C: ' num2str(mychans(c)) 'Mean P(' num2str(i(1)) '-' num2str(i(2)) ')'])
        end
    end

    text(.1,.1,['Mean P(' num2str(flim(1)) '-' num2str(flim(end)) ') given 10Hz Osc. Power dist(SU,HB) VS Vel(SL)'])
    ylabel('log10(P(10Hz)) of Spine to Head Distance')
    xlabel('Head Speed (cm/s)')


    A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@nanmean,nan);
    figure,imagescnan({vedges,sedges,clip(A,0,.015)'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ %tbp_phase.resample(yad);
% $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_median,nan);
% $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_var,nan);
% $$$ figure,imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,


% $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4|(fl>12&fl<18),chan),2);
% $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4,chan),2);



    %% End - Figure 6


    %% Figure 7 - rhm vs thetaM
    MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');
    %sname = 'jg05-20120317';
    sname = 'jg05-20120310';
    %sname = 'jg05-20120309';
    chans = [71:2:96];
    marker = 'spine_lower'

    Trial = MTATrial(sname,'all');
    Trial.ang.load(Trial);
    Trial.xyz.load(Trial);
    Trial.lfp.load(Trial,chans);

    wlfp = WhitenSignal(Trial.lfp.data,[],1);

    yl=[];
    for i = 1:Trial.lfp.size(2),
        [yl(:,:,i),fl,tl] = mtchglong(wlfp(:,i),2^12,Trial.lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
    end
    %yh=[];
    %for i = 1:Trial.lfp.size(2),
    %    [yh(:,:,i),fh,th] = mtchglong(wlfp(:,i),2^9,Trial.lfp.sampleRate,2^8,2^8*0.875,[],[],[],[40,120]);
    %end

    bang = ButFilter(Trial.ang(:,4,5,3),3,[1,20]./(Trial.ang.sampleRate./2),'bandpass');
    %bhh = ButFilter(Trial.xyz(:,7,3),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
    %bhx = ButFilter(Trial.xyz(:,7,1),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
    %bhy = ButFilter(Trial.xyz(:,7,2),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
    %if display, figure,plot([bang,bhh,bhx,bhy]+3.*repmat(1:4,size(bang,1),1)), end

    %wang = WhitenSignal([bang,bhh,bhx,bhy],[],1);
    wang = WhitenSignal([bang],[],1);
    [ya,fa,ta,phia,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^8,2^8*0.875,[],[],[],[2,16]);
    yad = MTADlfp([],[],ya,1/diff(ta(1:2)));
    yld = MTADlfp([],[],yl,1/diff(tl(1:2)));
    yld.resample(yad);
    xyz = Trial.xyz.copy;
    xyz.filter(gausswin(31)./sum(gausswin(31)));
    v = MTADxyz([],[],sqrt(sum(diff(xyz(:,marker,[1,2])).^2,3)).*Trial.xyz.sampleRate./10,Trial.xyz.sampleRate);
    v.resample(yad);

    count = 1;
    states = 'tvrwgl';
    for c = 1:13,
        for s = 1:numel(states)

            wper = Trial.stc{states(s)}.copy;
            wper.cast('TimeSeries');
            wper.resample(yad);



            ylpf = yld(:,fl<14&fl>5,c);yapf = yad(:,fa<14&fa>5,1,1);
            %yad(:,fa<14&fa>5,3,3)+yad(:,fa<14&fa>5,4,4)+yad(:,fa<14&fa>5,1,1);
            flpf = fl(fl<14&fl>5);fapf = fa(fa<14&fa>5);
            [ylpfsp,ylpfs] = sort(ylpf,2,'descend');[yapfsp,yapfs] = sort(yapf,2,'descend');
            flf = flpf(ylpfs(:,1));faf = fapf(yapfs(:,1));

            aind = wper.data==1;
            %aind = true(size(wper.data));
            lbounds = prctile(log10(ylpfsp(:,1)),[.5,99.5]);
            abounds = prctile(log10(yapfsp(:,1)),[.5,99.5]);
            aind = aind&~isinf(log10(ylpfsp(:,1)))&~isinf(log10(yapfsp(:,1)));
            %aind = aind&~isinf(log10(ylpfsp(:,1)))&~isinf(log10(yapfsp(:,1)))...
            %        &log10(ylpfsp(:,1))>lbounds(1)&log10(ylpfsp(:,1))<lbounds(2)...
            %        &log10(yapfsp(:,1))>abounds(1)&log10(yapfsp(:,1))<abounds(2);

            %figure,hist(flf,29)
            %figure,hist(faf,38)
            %figure,hist2([flf(aind),faf(aind)],29,38);caxis([0,140])

            %figure,hist2([clip(log10(ylpfsp(aind,1)),2.4,4), ...
            %              clip(log10(yapfsp(aind,1)),-5.5,-2.5)],30,30);caxis([0,130])

            %subplot(13,6,count),
            %figure,
            %plot(log10(ylpfsp(aind,1)),log10(yapfsp(aind,1)),'.')
            for shift = 1:50;
                [rho(c,s,shift),p(c,shift)] = corr(circshift(log10(ylpfsp(aind,1)),25-shift),log10(yapfsp(aind,1)));
            end
            %legend(num2str([rho,p]))
            %count = count+1;
        end
    end
    figure,imagesc(rho)


    figure,
    for i= 1:25
        subplot(5,5,i),hist2([circshift(flf,-i),faf],29,38);caxis([0,150])
    end

    %% End - Figure 7

  case 'fet_select'

    %% Figure 8 Feature Selection 
    MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');
    %sname = 'jg05-20120317';
    sname = 'jg05-20120310';
    %sname = 'jg05-20120309';
    display = 0;
    Trial = MTATrial(sname);
    Trial.load('xyz');
    Trial.load('ang');

    Trial.xyz.filter(gausswin(61)./sum(gausswin(61)));
    vel = Trial.vel;
    vel = clip(log10(vel),-2,2);
    vel =  MTADxyz('data',[zeros([1,size(vel,2)]);vel],'sampleRate',Trial.xyz.sampleRate);

    m=1;
    figure
    %hist2([vel(:,m),vel(:,7)],linspace(-2,2,64),linspace(-2,2,64));
    hist2([vel(vel(:,1)>0,m),vel(vel(:,1)>0,7)],linspace(-2,2,64),linspace(-2,2,64));
    caxis([0,600])


    % marker speed vs marker height
    figure
    ind = ':';
    vm = 1;
    hm = 1;
    hist2([vel(:,vm),clip(log10(Trial.xyz(:,hm,3)),0,3)],linspace(-2,2,64),linspace(1,2,64))
    caxis([0,1600])
    hold on
    states = 'rwgl';
    colors = 'wgmc';
    for i = 1:numel(states),
        ind = Trial.stc{states(i)};
        vh = [vel(ind,vm),clip(log10(Trial.xyz(ind,hm,3)),0,3)];
        eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
    end


    % marker speed vs marker segment pitch
    figure
    ind = ':';
    vm = 1;
    ms = [3,4];
    hist2([vel(:,vm),Trial.ang(:,ms(1),ms(2),2)],linspace(-1.5,2,64),linspace(-1.8,1.8,64))
    %ticks_lin2log
    caxis([0,1600])
    hold on
    states = 'rwgl';
    colors = 'wgmc';
    for i = 1:numel(states),
        ind = Trial.stc{states(i)};
        vh = [vel(ind,vm),Trial.ang(ind,ms(1),ms(2),2)];
        eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
    end



    %% End - Figure 8

    %% Figure - 9 Comodugram rhm and lfp

    %sname = 'jg05-20120317';
    %sname = 'jg05-20120315';
    sname = 'jg05-20120310.cof.all';
    %sname = 'jg05-20120309';
    %sname = 'jg04-20120129';
    %sname = 'jg04-20120130';
    Trial.lfp.filename = [Trial.name '.lfp'];
    chans = 68:3:95;

    Trial = MTATrial.validate(sname);

    %Trial.ang.load(Trial);
    stc = Trial.load('stc','MTAC_BATCH-fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-wrnpmsa');        
    xyz = Trial.load('xyz');
    lfp = Trial.load('lfp',chans);


    bang = fet_rhm(Trial,[],'mta');
    bang.resample(lfp);
    
    lbang = MTADlfp('data',WhitenSignal([lfp.data,bang.data],[],1),...
                    'sampleRate',lfp.sampleRate);

    states = {'gper','walk','rear','pause','groom','sit'};
    nsts = numel(states);
    nchan = numel(chans);
    figure,
    [Co,f] = Comodugram(Trial,lbang,stc.states(stc.gsi(states)));

    for s = 1:numel(states),
        for i =1:nchan
            subplot2(nchan,nsts,i,s);
            imagesc(f,f,Co(:,:,i,nchan+1,s)'),axis xy,
            if i==1,title(states{s}),end
            if i==1&s==1,ylabel([ 'Channel: ',num2str(chans(i))]),end
            caxis([-.5,.5])
            colormap jet
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        end
    end



  case 'rhm_distrb'

    %% Figure 11 - RHM distributions



    MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');

    %sname = 'jg05-20120315';
    %sname = 'jg05-20120309';
    %sname = 'jg04-20120129';
    %sname = 'jg04-20120130';
    %sname = 'co01-20140222';
    %sname = 'jg05-20120317';

    sname = 'jg05-20120310';
    Trial = MTATrial(sname,'all');
    %Trial.ang.load(Trial);
    Trial.xyz.load(Trial);
    Trial.xyz.filter(gausswin(9)./sum(gausswin(9)));

    rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
    hcom = Trial.com(rb);
    Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
    Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(gausswin(61)./sum(gausswin(61)),hcom),[1,3,2]));

    ang = Trial.ang.copy;
    ang.create(Trial);
    ang.data = [ang(:,4,5,3),ang(:,7,11,3)];
    %ang.resample(Trial.lfp);
    ang = ang.data;
    ang(isnan(ang(:,1)),:)=repmat(nanmean(ang),[sum(isnan(ang(:,1))),1]);

    wang = WhitenSignal(ang,[],1);

    [ya,fa,ta,pha,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^8, ...
                                    2^8*.875,[],[],[],[1,30]);



    chans = [71,81,86,96];
    Trial.lfp.load(Trial,chans);
    wlfp = WhitenSignal(Trial.lfp.data,[],1);
    yl=[];
    for i = 1:Trial.lfp.size(2),
        [yl(:,:,i),fl,tl] = mtchglong(wlfp(:,i),2^12,Trial.lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
    end



    figure,
    sp = [];
    for i = 1:numel(chans)
        sp(i) = subplot(6,1,i);
        imagesc(tl+2^11/Trial.lfp.sampleRate,fl,log10(yl(:,:,i)'));axis xy,caxis([1,3])
    end
    sp(6) = subplot(6,1,5);
    imagesc(ta+2^8/Trial.ang.sampleRate,fa,log10(ya(:,:,1,1)'));axis xy,caxis([-5,-2])
    sp(7) = subplot(6,1,6);
    imagesc(ta+2^8/Trial.ang.sampleRate,fa,log10(ya(:,:,2,2)'));axis xy,caxis([-5,-2])
    linkaxes(sp,'xy');



    figure,imagesc(tl,fl,log10(yl(:,:,1)'));axis xy




  case 'time_shifted_comodugram_rhm_lfp'
    %% Figure - 12 Time Shifted Comodugram rhm and lfp

    %MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');
    sname = 'jg05-20120310.cof.all';

    sname = 'jg05-20120317.cof.all';
    chans = 65:3:96;
    
    Trial = MTATrial.validate(sname);

    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);

    rhm = fet_rhm(Trial,[],'mta');
    rhm.resample(lfp);
    
    lbang = WhitenSignal([lfp.data,rhm.data],[],1);

    tshift = [-400,-250,-175,-70,-20,0,20,70,175,250,400];
    nts = numel(tshift);
    nchan = numel(chans);
    figure,
    sts = 'w';
    for s = 1:nts
        lb = MTADlfp([],[],[lbang(:,1:end-1),circshift(lbang(:,end),tshift(s))],lfp.sampleRate);

        lb = MTADlfp([],[],[lbang(:,3),circshift(lbang(:,end),tshift(s))],lfp.sampleRate);
        [Co,f] = Comodugram(Trial,lb,Trial.stc{sts});

        for i =1:nchan
            subplot2(nchan,nsts,i,s);
            imagesc(f,f,Co(:,:,i,nchan+1)');
            axis xy;
            if i==1,title(Trial.stc{sts}.label),end
            if i==1&s==1,ylabel([ 'Channel: ',num2str(chans(i))]),end
            caxis([-.65,.65])
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

        end

    end

    
    
  case 'state_dependent_comodugram_rhm_lfp'
    %% Figure - 12 Time Shifted Comodugram rhm and lfp

    %MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');
    sname = 'jg05-20120310.cof.all';

    sname = 'jg05-20120317.cof.all';
    chans = 65:3:96;
    
    Trial = MTATrial.validate(sname);

    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);

    rhm = fet_rhm(Trial,[],'mta');
    rhm.resample(lfp);
    
    lbang = WhitenSignal([lfp.data,rhm.data],[],1);




    figure,
    sts = 'awrpms';
    nsts = numel(sts);
    nchan = numel(chans);    
    for s = 1:nsts
        %lb = MTADlfp([],[],lbang(:,[1,end]),lfp.sampleRate);
        lb = MTADlfp([],[],lbang,lfp.sampleRate);
        [Co,f] = Comodugram(Trial,lb,Trial.stc{sts(s)});

        for i =1:nchan
            subplot2(nchan,nsts,i,s);
            imagesc(f,f,Co(:,:,i,nchan+1)');
            axis xy;
            if i==1,title(Trial.stc{sts(s)}.label),end
            if i==1&s==1,ylabel([ 'Channel: ',num2str(chans(i))]),end
            caxis([-.65,.65])
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

        end

    end
    

  case 'marker_reconstruction noise'
    %% Figure 13 - Marker to marker distances
    MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    %MTAConfiguration('/data/data/gravio','absolute');
    sname = 'jg05-20120310';
    Trial = MTATrial(sname,'all');
    Trial.load('ang');
    Trial.load('xyz');

    figure
    for m = 1:8,
        for o = 1:8,
            subplot2(8,8,m,o)
            hist(Trial.ang(Trial.ang(:,m,o,3)~=0,m,o,3),500)
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        end
    end


  case 'bhv_lfp_psd'
    %% Figure 14 - Spectral power for each behavior
    %sname = 'jg05-20120317';
    sname = 'jg05-20120310';
    %sname = 'jg05-20120309';
    display = 0;
    %chans = [1:2:8];
    chans = [65:1:96];
    marker = 'spine_lower'

    Trial = MTATrial(sname,'all');
    Trial.ang.load(Trial);
    Trial.xyz.load(Trial);
    Trial.lfp.load(Trial,chans);
    Trial.stc.updateMode('auto_wbhr');
    Trial.stc.load;

    wlfp = WhitenSignal(Trial.lfp.data,[],1);

    matlabpool open 12

    tl=[];
    fl=[];
    yl=[];
    spectral.nfft = 2^11;
    spectral.window = 2^10;
    parfor i = 1:Trial.lfp.size(2),
        [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,Trial.lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],[1,40]);
    end
    fl = fl(:,1);
    tl = tl(:,1);
    yld = MTADlfp('data',yl,'sampleRate',1/diff(tl(1:2,1)));
    yld.data(yld.data==0)=nan;
    yld.data = log10(yld.data);
    yld.data = (yld.data-repmat(nanmedian(yld.data),[yld.size(1),1,1]))./repmat(nanstd(yld.data),[yld.size(1),1,1]);
    tshift = round(spectral.window/2/Trial.lfp.sampleRate*yld.sampleRate);




    sts='r'
    evt=1;
    figure,cnt=1;
    for i=linspace(-2,2,40).*yld.sampleRate
        imagesc(fl,1:13,sq(nanmean(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
        caxis([-1,1])
        text( 35,13,num2str((i)./yld.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
        pause(.2)
        frms(cnt) = getframe;
        cnt=cnt+1;
    end

    prm.fps = 5;
    prm.loop = inf;
    makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_on_theta.gif',prm);


    th=[];
    fh=[];
    yh=[];
    spectral.nfft = 2^9;
    spectral.window = 2^7;
    spectral.freq = [40,120];
    parfor i = 1:Trial.lfp.size(2),
        [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,Trial.lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],spectral.freq);
    end
    fh = fh(:,1);
    th = th(:,1);
    yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2,1)));
    yhd.data(yhd.data==0)=nan;
    yhd.data = log10(yhd.data);
    yhd.data = (yhd.data-repmat(nanmedian(yhd.data),[yhd.size(1),1,1]))./repmat(nanstd(yhd.data),[yhd.size(1),1,1]);
    tshift = round(spectral.window/2/Trial.lfp.sampleRate*yhd.sampleRate);

    sts='r'
    evt=2;
    figure,cnt=1;
    for i=linspace(-2,2,200).*yhd.sampleRate
        imagesc(fh,1:yhd.size(3),sq(nanmean(yhd(round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
        caxis([-.51,.51])
        text( 80,30,num2str(i./yhd.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
        pause(.2)
        frms(cnt) = getframe;
        cnt=cnt+1;
    end

    prm.fps = 10;
    prm.loop = inf;
    makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_off_gamma.gif',prm);


    figure,
    s = 'r'
    sp1 = subplot(121);hold on
    sp2 = subplot(122);hold on
    colors= jet(numel(chans))';

    for c = colors
        plot(sp1,fl,-.5*find(ismember(colors',c','rows'))+nanmean(log10(yld(Trial.stc{s,yld.sampleRate}.data(2:end-1,:),:,ismember(colors',c','rows')))),'color',c)
        plot(sp2,fh,-.21*find(ismember(colors',c','rows'))+nanmean(log10(yhd(Trial.stc{s,yhd.sampleRate},:,ismember(colors',c','rows')))),'color',c)
    end


    figure,
    %subplot(121),
    chan = 15;
    st1 = 'l'; st2 = 'g'; 
    st1 = 'r'; st2 = 'w';
    boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-b','alpha')

    figure,
    %subplot(121),
    st1 = 'l'; st2 = 'g';
    boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-b','alpha')


    subplot(122),

    boundedline(fh,mean(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),'-r',fh,mean(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),'-b','alpha')


    evt = 1;
    sts = 'r';
    figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<12&fl>6,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')

    evt = 2;
    sts = 'r';
    figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<25&fl>18,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')


    evt = 1;
    sts = 'r';
    figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yhd.data(:,fh<120&fh>80,:),2)),round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yhd.sampleRate),round(4*yhd.sampleRate),0),2))')

    tdr = sq(mean(yl(:,fl>6&fl<12,:),2)./(mean(yl(:,fl<5,:),2)+mean(yl(:,fl<18&fl>13,:),2))/2);
    tdr = MTADlfp('data',tdr,'sampleRate',1/diff(tl(1:2,1)));


    tdr.data = (tdr.data-repmat(nanmedian(tdr.data),[tdr.size(1),1,1]))./repmat(nanstd(tdr.data),[tdr.size(1),1,1]);


    evt = p1;
    sts = 'r';
    figure,imagesc(sq(nanmean(GetSegs(tdr.data,round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate-tshift/yld.sampleRate),round(4*yld.sampleRate),0),2))')

    sts = 'w';
    [U,V,D] = svd(cov(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-tshift/yld.sampleRate),:,21)));
    figure,imagesc(fl,fl,U),axis xy



  case 'uhhh'
    %% Figure # Unit state place field formation

    sname = 'jg05-20120317.cof.all';
    Trial = MTATrial.validate(sname);
    Trial.stc.updateMode('auto_wbhr');
    Trial.stc.load;
    states = 'rwgl';

    for i = 1:numel(states),
        spk{i} = Trial.spk.copy;
        spk{i}.create(Trial,Trial.xyz.sampleRate,states(i));        
        sts{i} = Trial.stc{states(i)}.copy;
        sts{i}.cast('TimeSeries');
        endn

        clrs = 'rbcm';

        unit = 71;
        figure
        hold on
        for i = 1:numel(states),
            Lines(spk{i}(unit),[],clrs(i));
        end
        for i = 1:numel(states),
            plot(sts{i}.data*2+i/4,clrs(i))
        end
        ylim([-3,4])

        plot([sqrt(sum(sq([Trial.xyz(:,7,1)-75,Trial.xyz(:,7,2)+50]).^2,2))<100]-2)

        %plot((Trial.xyz(:,7,3)-Trial.xyz(:,1,3))./max(Trial.xyz(:,7,3)-Trial.xyz(:,1,3)),'r')

    end

    
    
    
    
% case 'head_motion' -----------------------------------------------------------------------------
  case 'head_motion'
    Trial = MTATrial('jg05-20120317');
    fs = []; ts = [];

    xyz = Trial.load('xyz');
    % create a ridgid body model
    rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
    % find the center of mass of the model
    hcom = xyz.com(rb);
    % add coordinates of the model's center of mass to the xyz object
    xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));

    xyz.filter('ButFilter',3,50,'low');

    ang = create(Trial.ang.copy,Trial,xyz);
    bang = ButFilter(ang(:,'head_back','fhcom',3),3,[2,30]./(Trial.ang.sampleRate/2),'bandpass');
    bang = [bang,ButFilter(ang(:,'head_right','fhcom',3),3,[2,30]./(Trial.ang.sampleRate/2),'bandpass')];
    bang = [bang,ButFilter(ang(:,'head_top','fhcom',3),3,[2,30]./(Trial.ang.sampleRate/2),'bandpass')];
    bang = [bang,circ_dist(ang(:,'head_back','head_front',2),ang(:,'head_back','fhcom',2))];
    bang = WhitenSignal(bang,[],1);
    
    figure,plot(circ_dist(ang(:,'head_back','head_front',2),ang(:,'head_back','fhcom',2)))

    bfet = Trial.xyz.copy;
    bfet.data = bang;

    [ys,fs,ts,phi,fst] = fet_spec(Trial,bfet,'mtchglong',true,'overwrite',true);

    xyz.resample(ys);
    wang = create(Trial.ang.copy,Trial,xyz);    
    
    c = 4;
    f = 35;
    nind = nniz(wang(:,1,2,1));
    figure,hist2([wang(nind,5,7,2),log10(ys(nind,f,c,c))],linspace(-1.5,1.5,100),linspace(-7,-2,100));
    s = 'w';
    nind = Trial.stc{s};
    figure,hist2([wang(nind,5,7,2),log10(ys(nind,f,c,c))],linspace(-1.5,1.5,100),linspace(-7,-2,100));
    
    
    nc = ys.size(3);
    figure,sp = [];
    for i = 1:nc,
        sp(i) = subplot(nc,1,i);imagesc(ts,fs,log10(ys(:,:,i,i))'),caxis([-7,-2]),axis xy
    end
    linkaxes(sp,'xy');

    c = 2;
    f = 35;
    figure,hist2(log10([ys(:,35,1,1),ys(:,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'m'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'w'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'r'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    
    vfet = vel(Trial.load('xyz')
    [ys,fs,ts,phi,fst] = fet_spec(Trial,bfet,'mtchglong',true,'overwrite',true);
    
     pbins = linspace(-7,-3,100);
     pfd = histc(log10(ys(:,:,1,1)),pbins,1);
     figure,imagesc(pbins,fs,pfd'),axis xy


     hfig = figure(848283); 
     nind = nniz(ys);
     for s = 'arwms';
         clf(hfig)
         nind = Trial.stc{'w'};
         pfd = histc(log10(ys(nind,:,1,1)),pbins,1);
         imagesc(pbins,fs,pfd'),axis xy
         saveas(
     end
     
     
     fet =e fet_lgr(Trial);
     nind = nniz(fet);

     
     
     
     
     [isig] = fastica(cov(fet(nind,:)));
     nind = nniz(fet);figure,N = hist2([fet(nind,:)*isig(11,:)',fet(nind,:)*isig(2,:)'],linspace(-10,10,100),linspace(-15,15,100));
     


     
% case 'dz-lfp comodulation' ------------------------------------------------------------------------
  case 'dz-lfp comodulation'
    
    %% comodulation of low frequency oscillations in the lowerspine
    %% marker and the lfp
    Trial = MTATrial('jg05-20120317');
    Trial.load('stc','hand_labeled_rev2');
    xyz = Trial.load('xyz');
    fxyz = xyz.copy;
    fxyz.filter('ButFilter',3,40,'low');
    ang = create(MTADang,Trial,fxyz);

    lfp = Trial.load('lfp',84);
    name = 'lower spine Z speed'; label = 'lszs'; key = 'z';
    zv = MTADfet.encapsulate(Trial,...
                         diff(ang(:,1,4,3)),...
                         xyz.sampleRate,...
                         name,label,key);

    zv = MTADfet.encapsulate(Trial,...
                         diff(circ_dist(ang(:,2,4,1),ang(:,4,7,1))),...
                         xyz.sampleRate,...
                         name,label,key);
%lfp.resample(zv);
    zv.resample(lfp);
    zv.filter('ButFilter',3,40,'low');


    dspec = struct('nFFT',      2^10,...
                   'SampleRate',zv.sampleRate,...
                   'WinLength', 2^9,...
                   'FreqRange', [1,150]);
    dspec = struct2varargin(dspec);
    nind = nniz(lfp)&nniz(zv);
    nind = Trial.stc{'w+n+p'};
    [Co,f] = Comodugram([lfp(nind),zv(nind)],dspec{:});

    figure;
    for xi = 1:2,
        for yi = 1:2,
            subplot2(2,2,xi,yi);
            imagesc(f,f,Co(:,:,xi,yi)');
            axis xy;
            colormap jet
            caxis([0,1])
        end
    end
    
    zvp = zv.copy;
    zvp = zvp.phase([2,4]);
    spk = Trial.spk.create(Trial,Trial.xyz.sampleRate,'walk');

    thp = lfp.phase([6,14]);
    figure,rose(zvp(spk(10)),30)

    figure
    for sid = 1:90,
        clf;
        plot(thp(spk(sid)),zvp(spk(sid)),'.');
        title(num2str(sid));
        pause(.4)
    end
    
% case 'spk_trig_lfp_ave'    ---------------------------------------------------------------
  case 'spk_trig_lfp_ave'
    sname = 'jg05-20120310.cof.all';

    sname = 'jg05-20120317.cof.all';
    chans = 65:1:96;
    
    Trial = MTATrial.validate(sname);

    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);
    
    sts = 'awrpms';
    nsts = numel(sts)
    for s = 1:nsts,
        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,xyz.sampleRate,sts(s),[],'deburst');
        for u = spk{s}.map(:,1)',
            try,
                spkTrigAve(:,:,u,s) = sq(mean(lfp.segs(spk{s}(u)-120,240),2));
            end
        end
    end

    hfig = figure,
    for u = 1:size(spkTrigAve,3)
        clf
        for s = 1:nsts,
            subplot(nsts,1,s);
            imagesc(spkTrigAve(:,:,u,s)');
            title(['unit: ' num2str(u) sts(s)]);
            caxis([-3000,3000]);
            colormap jet;
        end
        reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
          hfig,                                ... Figure handle
          ['lfpSpkTrigAve'],                   ... Figure Set Name
          'req',                               ... Directory where figures reside
          false,                               ... Do Not Preview
          [Trial.filebase '-' num2str(u)],     ... Tumbnail caption
          [Trial.filebase '-' num2str(u)],     ... Expanded caption
          [],                                  ... Resolution
          false,                               ... Do Not Save FIG
          'png',4,8);                                % Output Format
    end

  case 'spk_trig_lfp_psd_ave'
    sname = 'jg05-20120310.cof.all';
    sname = 'jg05-20120317.cof.all';

    chans = 65:1:96;
    
    Trial = MTATrial.validate(sname);

    Trial.stc.load(Trial,['MTAC_BATCH-fet_mis'...
                          '_SR_12_NORM_1_REF_jg05-20120317.cof.all'...
                          '_STC_hand_labeled_rev3_jg_NN_100_NI_100'...
                          '_NN_multiPN_RAND_WSBNT-wrnpmsa']);

    xyz = Trial.load('xyz');
    
    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);
    lfp.resample(xyz);
    
    sts = 'awrpms';
    nsts = numel(sts)
    for s = 6:nsts,
        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,Trial.xyz.sampleRate,sts(s),[],'deburst');
        for u = spk{s}.map(:,1)',
            try,
                spkTrigAve(:,:,u,s) = sq(mean(lfp.segs(spk{s}(u)-60,120),2));
            end
        end
    end

    hfig = figure,
    for u = 1:size(spkTrigAve,3)
        clf
        for s = 1:nsts,
            subplot(nsts,1,s);
            imagesc(spkTrigAve(:,:,u,s)');
            title(['unit: ' num2str(u) sts(s)]);
            caxis([-3000,3000]);
            colormap jet;
        end
        reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
          hfig,                                ... Figure handle
          ['lfpSpkTrigAve'],                   ... Figure Set Name
          'req',                               ... Directory where figures reside
          false,                               ... Do Not Preview
          [Trial.filebase '-' num2str(u)],     ... Tumbnail caption
          [Trial.filebase '-' num2str(u)],     ... Expanded caption
          [],                                  ... Resolution
          false,                               ... Do Not Save FIG
          'png',4,8);                                % Output Format
    end


    
    
% case 'spk_trig_lfpPSD_ave' -----------------------------------------------------------------
  case 'spk_trig_lfpPSD_ave'
    sname = 'jg05-20120310.cof.all';
    sname = 'jg05-20120317.cof.all';

    chans = 65:1:96;
    
    Trial = MTATrial.validate(sname);

    Trial.stc.load(Trial,['MTAC_BATCH-fet_mis'...
                          '_SR_12_NORM_1_REF_jg05-20120317.cof.all'...
                          '_STC_hand_labeled_rev3_jg_NN_100_NI_100'...
                          '_NN_multiPN_RAND_WSBNT-wrnpmsa']);

    xyz = Trial.load('xyz');
    
    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);

    wlfp = lfp.copy;
    wlfp.data = WhitenSignal(lfp.data,[],1);
    
    
    parspec = struct('nFFT',2^10,...
                     'Fs',  wlfp.sampleRate,...
                     'WinLength',2^9,...
                     'nOverlap',2^9*.5,...
                     'NW',3,...
                     'Detrend',[],...
                     'nTapers',[],...
                     'FreqRange',[1,140]);
    

    swlfp = wlfp.copy;
    swlfp.data = wlfp(:,1);
    [ys,fs,ts] = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);   
    for c = 2:numel(chans),
        swlfp.data = wlfp(:,c);
        tys = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);
        ys.data(:,:,c) = tys.data;
    end

        
    nys = ys.copy;
    nys.unity;

    sts = 'awrpms';
    nsts = numel(sts)
    for s = 1:nsts,
        spk{s} = Trial.spk.copy;
        spk{s}.create(Trial,ys.sampleRate,sts(s),[],'deburst');
    end
    
    clear('spkTrigAve');
    for s = 1:nsts,
        for u = spk{s}.map(:,1)',
            try,
                spks = spk{s}(u);
                spkTrigAve(:,:,u,s) = sq(nanmean(nys(spks,:,:)));
            end
        end
    end

    hfig = figure(30239230);
    for u = 1:size(spkTrigAve,3)
        clf
        for s = 1:nsts,
            subplot(nsts,1,s);
            imagesc(fs,chans,spkTrigAve(:,:,u,s)');
            title(['unit: ' num2str(u) sts(s)]);
            caxis([-1,1]);
            colormap jet;
        end
        reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
          hfig,                                ... Figure handle
          ['lfpPSDSpkTrigAve'],                   ... Figure Set Name
          'req',                               ... Directory where figures reside
          false,                               ... Do Not Preview
          [Trial.filebase '-' num2str(u)],     ... Tumbnail caption
          [Trial.filebase '-' num2str(u)],     ... Expanded caption
          [],                                  ... Resolution
          false,                               ... Do Not Save FIG
          'png',3,9);                                % Output Format
    end

    
% case 'get_back_on_task'    --------------------------------------------------------------
    case 'get_back_on_task'
    
    sname = 'jg05-20120310.cof.all';
    sname = 'jg05-20120317.cof.all';

    chans = 65:1:96;
    
    Trial = MTATrial.validate(sname);

    Trial.stc.load(Trial,['MTAC_BATCH-fet_mis'...
                          '_SR_12_NORM_1_REF_jg05-20120317.cof.all'...
                          '_STC_hand_labeled_rev3_jg_NN_100_NI_100'...
                          '_NN_multiPN_RAND_WSBNT-wrnpmsa']);

    xyz = Trial.load('xyz');
    
    Trial.lfp.filename = [Trial.name,'.lfp']; % Remove later
    lfp = Trial.load('lfp',chans);

    wlfp = lfp.copy;
    [~,arm] = WhitenSignal(lfp(nniz(lfp.data),1));
    wlfp.data = WhitenSignal(lfp.data,[],1,arm);    
    
    n = 10;    
    frng = fliplr(2:1:145);

    f = frng(1);    

    wb = round(1/f*lfp.sampleRate*n);        
    parspec = struct('nFFT',wb*2,...
                     'Fs',  wlfp.sampleRate,...
                     'WinLength',wb,...
                     'nOverlap',round(wb*.875),...
                     'NW',3,...
                     'Detrend',[],...
                     'nTapers',[],...
                     'FreqRange',[f-1,f+1]);
    

    swlfp = wlfp.copy;
    swlfp.data = wlfp(:,1);
    [ys,fs,ts] = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);    

    fsa = [];
    fsa = cat(1,fs,fsa);        


    for f = frng(2:end),
        wb = round(1/f*lfp.sampleRate*n);        
        parspec = struct('nFFT',wb*2,...
                         'Fs',  wlfp.sampleRate,...
                         'WinLength',wb,...
                         'nOverlap',round(wb*.875),...
                         'NW',3,...
                         'Detrend',[],...
                         'nTapers',[],...
                         'FreqRange',[f-5,f+5]);
        [tys,fs] = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);
        fsa = cat(1,fs,fsa);        
        tys.resample(ys);
        ys.data = cat(2,tys.data,ys.data);
    end

    ysg = ys.copy;

    [fsg,fid] = sort(fsa);
    ysg.data = ysg(:,fid);
    figure,imagesc(ts,1:846,log10(abs(ysg.data))')    
    caxis([0,4]),colormap jet    

    figure,
    imagesc(ts,1:846,nunity(log10(abs(ysg.data)))')    
    caxis([-3,3]),colormap jet    
    axis xy
    
    [fsg,fid] = sort(fsa);    
    [fsu,fud] = unique(fsg);

    ysn = ysg.copy;
    ysn.data = ysn(1:2000,fud);    
    ysn.data = interp1(ysn.data',repmat(fsu,1,ysn.size(1))',repmat([.5:.5:140],ysn.size(1),1)','spline');    

    ysp = ysn.copy;
    ysp.data = zeros([size(ysn,1),numel([.5:.5:140])]);    
    for t = 1:size(ysn,1),
        try
        ysp.data(t,:) = interp1(ysn(t,:),fsu',[.5:.5:140],'linear'); 
        end
    end
    
    f = 140
    swlfp = wlfp.copy;
    swlfp.data = wlfp(:,1);
    [ys,fs,ts] = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);   
    for c = 2:numel(chans),
        swlfp.data = wlfp(:,c);
        tys = fet_spec(Trial,swlfp,'mtchglong',false,[],parspec);
        ys.data(:,:,c) = tys.data;
    end
    
  case 'ufr copulas'
    Trial = MTATrial.validate('jg05-20120310');
    
    xyz = Trial.load('xyz');
    ufr = create(MTADufr,Trial,xyz,[],[],0.2);

    u = [29,24];
    figure,
    hist2([MakeUniformDistr(ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(1))),...
           MakeUniformDistr(ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(2)))],...
          20,...
          20);

    figure
    u = [11,32];
    hist2(log10(abs([ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(1))+randn(size(ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(1)),1),1),...
           ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(2))+randn(size(ufr(ufr(:,u(1))~=0|ufr(:,u(2))~=0,u(1)),1),1)])),...
          linspace(1,2,100),...
          linspace(1,2,100));
    
  case 'BHVPFS'
    Trial= MTATrial('jg05-20120310');    
    Trial.load('stc','nn0317_PP');
    states = Trial.stc.list_state_attrib;
    binDims = [40,40];
    smoothingWeights = [1.2,1.2];
    units = [];
    overwrite = false;
    for s = 1:numel(states)
        pfs{s} = MTAApfs(Trial,units,states{s},overwrite, ...
                         'binDims',binDims,'SmoothingWeights',smoothingWeights);
    end
    units = pfs{1}.data.clu;    

    [accg,tbin] = autoccg(Trial,units,'theta');

    t = 1;
    mRate = pfs{9}.maxRate;
    %units = find(sq(mRate(1,1,:))>3);



    autoincr = false;
    hfig = figure(849274);    
    unit = units(1);
    while unit~=-1,
        for s = 1:numel(states)
            subplot(3,4,s)
            hold('on')
            pf = pfs{s};
            ratemap = pf.plot(unit,'isCircular',true);
            ratemap(isnan(ratemap)) = -1;
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
            text(pf.adata.bins{1}(1)+30,pf.adata.bins{2}(end)-50,...
                 sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
            colormap([0,0,0;parula]);


% $$$             plot(peakPatchCOM(t,i,unit==units,1),...
% $$$                  peakPatchCOM(t,i,unit==units,2),'*k');
% $$$             xlim([-600,600]),ylim([-350,350])                    
            title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);
        end   
        ForAllSubplots('colormap([0,0,0;parula])');
        ForAllSubplots(['caxis([-1,',num2str(sq(mRate(unit))),'.*1.5])']);

        subplot(3,4,10)
        bar(tbin,accg(:,unit));
        xlim([min(tbin),max(tbin)]);
        title([' AutoCCG: Unit ',num2str(unit)]);

        unit = figure_controls(hfig,unit,units,autoincr);
    end
    
    
    
    
    
    
end











