sli = 1
rli = 1;

rlist = SessionList('training_hand_labeled');
slist = {'hand_labeled_jg';'hand_labeled_Ed'};

for rli = 1:numel(rlist),

    
    fetSet  = 'fet_tsne_rev13';
    sampleRate = 12;
    nNeurons = 100;
    nIter = 200;
    states = {'walk','rear','turn','pause','groom','sit'};
    rndMethod = 'WSBNT';
    norm = true;
    mref = true;



    model = ['MTAC_BATCH-' fetSet ...
             '_SR_' num2str(sampleRate) ...
             '_NORM_' num2str(norm) ...         
             '_REF_' rlist(rli).sessionName, '.' rlist(rli).mazeName '.' rlist(rli).trialName ...
             '_STC_' rlist(rli).stcMode ...
             '_NN_' num2str(nNeurons) ...
             '_NI_' num2str(nIter) ...         
             '_NN_multiPN_RAND_' rndMethod];
    
    
    refSession = MTATrial.validate(rlist(rli));
    rfet = feval(fetSet,refSession,sampleRate,false);
    [~,refMean,refStd] = nunity(rfet(refSession.stc{'a'},:));


    for sli = 1:numel(slist),

        SesList = SessionList(slist{sli});
        mapped = '-map2ref';
        ds = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,mapped,'.mat']));

        %% Pre-process features
        for s = 1:numel(SesList);
            Trial = MTATrial.validate(SesList(s));

            features = feval(fetSet,Trial,sampleRate,false);
            features.map_to_reference_session(Trial,refSession);
            features.unity([],refMean,refStd);

            if isempty(Trial.fet),
                Trial.fet = MTADfet(Trial.spath,...
                                    [],...
                                    [],...
                                    [],...
                                    Trial.sync.copy,...
                                    Trial.sync.data(1),...
                                    []);                  
            end


            xyz = Trial.load('xyz');

            try 
                ss = Trial.load('fet','3dss');
                ss.resample(xyz);
            catch
                txyz = xyz.data;
                pnts = zeros([xyz.size(1),105,3]);
                for ind = 1:xyz.size(1),
                    try
                        pnts(ind,:,:) = fnplt(cscvn(sq(txyz(ind,1:4,:))'))';
                    end
                end
                name = '3d spline interpolated spine'; label = '3dss'; key = 's';
                ss = MTADfet.encapsulate(Trial,...
                                         pnts,...
                                         xyz.sampleRate,...
                                         name,label,key);
                ss.updateFilename(Trial);
                ss.save;
            end
            xyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);

            xyz.addMarker('bcom',[.7,0,.7],{},...
                          xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
            xyz.addMarker('hcom',[.7,0,.7],{},...
                          xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

            StcHL = Trial.stc.copy;
            shl = MTADxyz('data',double(~~stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
% $$$ 
% $$$ % $$$ cstc = 1-max(ds.d_state{s},[],2)./nIter;
% $$$ % $$$ figure,plot(cstc)
% $$$ % $$$ Lines([],0.05,'k');
% $$$ % $$$ Lines(StcHL{'w'}(:),[],'c');
% $$$ 
            ysm = MTADxyz('data',double(0<stc2mat(ds.stc{s},xyz,states)),'sampleRate',xyz.sampleRate); 
% $$$ 
% $$$ %ysm = MTADxyz('data',double(bsxfun(@rdivide,p_state{s},max(p_state{s},[],2))==1),'sampleRate',xyz.sampleRate); 
% $$$ 
% $$$ % $$$ figure,sp=[];
% $$$ % $$$ sp(end+1)=subplot(211);
% $$$ % $$$ imagesc(shl.data')
% $$$ % $$$ sp(end+1)=subplot(212);
% $$$ % $$$ imagesc(ysm.data')
% $$$ % $$$ linkaxes(sp,'xy')
% $$$ 
% $$$ labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');
% $$$ %ind = cstc<0.05;
% $$$ ind = any(ysm.data,2)&any(shl.data,2);
% $$$ 
% $$$ 
% $$$ tcm = confmat(shl(ind&labelingEpochs.data,:),ysm(ind&labelingEpochs.data,:)); % #DEP: netlab
% $$$ confusionMatrix = round(tcm./xyz.sampleRate,2);
% $$$ precision = round(diag(tcm)./sum(tcm,2),4).*100;
% $$$ sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
% $$$ accuracy = sum(diag(tcm))/sum(tcm(:));
% $$$ 
% $$$ precision'
% $$$ sensitivity
% $$$ accuracy

            ang = create(MTADang,Trial,xyz);

            %% COMPOSITE: rear
            StcCor = ds.stc{s}.copy;
            rhh = [];
            rthresh = -0.2;
            tds = ds.d_state{s};
            [tpv,tps] = sort(ds.p_state{s},2,'descend');
            try
            for rp = ds.stc{s}{'r'}.data',
                rhh(end+1) = min(((xyz(rp',10,3)-refMean(4))./refStd(4)).*((ang(rp',3,4,2)-refMean(8))./refStd(8)));
                %rhh(end+1) = max(xyz(rp',10,3).*ang(rp',3,4,2));
                if rhh(end)<rthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    StcCor.states{StcCor.gsi(states{mode(tps(pind,1))})}.data = ...
                        [StcCor.states{StcCor.gsi(states{mode(tps(pind,1))})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('r')}.data(rhh<rthresh,:) = [];
            end

            %% DISTANCE: Walk
            try
            wd = [];
            wthresh = 1.5;
            for rp = StcCor{'w'}.data',
                wd(end+1) = log10(sqrt(sum([xyz(rp(2),1,[1,2])-xyz(rp(1),1,[1,2])].^2,3)));
                if wd(end)<wthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    if msts == 2; msts = 4; end
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('w')}.data(wd<wthresh,:) = [];
            end

            
            %% SPEED: Pause
% $$$             figure,
% $$$             eds = linspace(-.2,1,40);
% $$$             pd = [];
% $$$             %for rp = StcHL{'s'}.data',            
% $$$             for rp = ds.stc{s}{'s'}.data',                            
% $$$                 pd(end+1) = circ_mean(ang(rp',2,3,2));
% $$$                 %pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
% $$$             end
% $$$             subplot(211),bar(eds,hist(pd,eds),'histc')
% $$$             pd = [];
% $$$             %for rp = StcHL{'p'}.data',            
% $$$             for rp = ds.stc{s}{'p'}.data',                                            
% $$$                 pd(end+1) = circ_mean(ang(rp',2,3,2));                
% $$$                 %pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
% $$$             end
% $$$             subplot(212),bar(eds,hist(pd,eds),'histc')


            %% HEIGHT: Pause
            try
            key = 'p';
            pd = [];
            pthresh = -0.37;
            for rp = StcCor{key}.data',
                pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
                if pd(end)<pthresh
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    msts = 6;
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
p                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
            end


            
            %% HEIGHT: Pause
            try
            key = 's';
            pd = [];
            pthresh = -0.37;
            for rp = StcCor{key}.data',
                pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
                if pd(end)>pthresh
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    msts = 6;
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd>pthresh,:) = [];
            end


            %% SPEED: pause
            try
            pd = [];
            pthresh = -0.1;
            key = 'p';
            %for rp = StcCor{'s'}.data',            
            for rp = StcCor{key,features.sampleRate}.data',
                %pd(end+1) = mean((xyz(rp',4,3)-refMean(2))./refStd(2));
                pd(end+1) = min(features(rp',5));    
                if pd(end)>pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = 4;%mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd>pthresh,:) = [];
            end


            

            %% DURATION: Pause
            try
            pd = [];
            pthresh = 1.3;
            key = 'p';
            for rp = StcCor{key}.data',
                pd(end+1) = log10(abs(diff(rp)));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    if msts == 2; 
                        tps(pind,:) = circshift(tps(pind,:),-1,2);
                        msts = mode(tps(pind,1));
                    end                    
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
            end
            
            

            %% Correct sit Duration
            %% DURATION
            try
            pd = [];
            pthresh = 2.4;
            for rp = StcCor{'s'}.data',
                %for rp = StcHL{'s'}.data',    
                pd(end+1) = log10(abs(diff(rp)));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = 4;%mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('s')}.data(pd<pthresh,:) = [];
            end




            %% DURATION: Turn
            try
            pd = [];
            ad = [];
            pthresh = 1.8;
            athresh = 0.3;
            key = 'n';
            for rp = StcCor{key}.data',
                pd(end+1) = log10(abs(diff(rp)));
                ad(end+1) = abs(circ_dist(ang(rp(2),1,4,1),ang(rp(1),1,4,1)));
                if pd(end)<pthresh&ad(end)<athresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    if msts == 2; 
                        tps(pind,:) = circshift(tps(pind,:),-1,2);
                        msts = mode(tps(pind,1));
                    end                    
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh&ad<athresh,:) = [];
            end


            
             
            
            %walk 
            try
            pd = [];
            ad = [];
            pthresh = 1.5;
            athresh = 0.3;
            key = 'w';
            for rp = StcCor{key}.data',
                %for rp = StcHL{key}.data',
                pd(end+1) = log10(abs(diff(rp)));
                ad(end+1) = abs(circ_dist(ang(rp(2),1,4,1),ang(rp(1),1,4,1)));
                if pd(end)<pthresh&ad(end)>athresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    if msts == 2; 
                        tps(pind,:) = circshift(tps(pind,:),-1,2);
                        msts = mode(tps(pind,1));
                    end                    
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh&ad>athresh,:) = [];
            end
            
            
            %% DURATION: Sit
            try
            pd = [];
            pthresh = 2.5;
            for rp = StcCor{'s'}.data',
                %for rp = StcHL{'s'}.data',    
                pd(end+1) = log10(abs(diff(rp)));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = 4;%mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('s')}.data(pd<pthresh,:) = [];
            end

            
            
            
            csm = [];
            csm = MTADxyz('data',double(0<stc2mat(StcCor,xyz,states)),'sampleRate',xyz.sampleRate); 
            labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

            errorEpochs = Trial.stc{'e'}.cast('TimeSeries');
            if isempty(errorEpochs)
                errorEpochs = labelingEpochs.copy;
                errorEpochs.data = ~errorEpochs.data;
            end
            
            %ind = cstc<0.05;
            ind = MTADepoch('data',any(csm.data,2)&any(shl.data,2),...
                            'sampleRate',xyz.sampleRate,...
                            'type','TimeSeries');
            ind.resample(labelingEpochs);
            errorEpochs.resample(labelingEpochs);
            ind = ind.data&labelingEpochs.data&~errorEpochs.data;

            tcm = confmat(shl(ind,:),csm(ind,:)); % #DEP: netlab
            ls{s}.confusionMatrix = round(tcm./xyz.sampleRate,2);
            ls{s}.precision = round(diag(tcm)./sum(tcm,2),4)'.*100;
            ls{s}.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
            ls{s}.accuracy = sum(diag(tcm))/sum(tcm(:));

            StcCor.updateMode([StcCor.mode,'_PP']);
            stc{s} = StcCor.copy;
        end



        mapped = '-map2ref';

        save(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'_PP',mapped,'.mat']),...
             '-v7.3','slist','rlist','nNeurons','nIter','sampleRate','model','fetSet','rndMethod',...
             'states','stc','d_state','p_state','ls');


    end
end




csm = [];
csm = MTADxyz('data',double(0<stc2mat(StcCor,xyz,states)),'sampleRate',xyz.sampleRate); 
figure,sp=[];
sp(end+1)=subplot(311);imagesc(shl.data')
sp(end+1)=subplot(312);imagesc(ysm.data')
sp(end+1)=subplot(313);imagesc(csm.data')
linkaxes(sp,'xy')


figure,sp=[];
sp(end+1)=subplot(511);imagesc(ds.p_state{s}');caxis([-100,100]);
sp(end+1)=subplot(512);imagesc(ds.d_state{s}');
sp(end+1)=subplot(513);imagesc(shl.data')
sp(end+1)=subplot(514);imagesc(ysm.data')
sp(end+1)=subplot(515);imagesc(csm.data')
linkaxes(sp,'xy')



for sli = 1:numel(slist)
    SesList = SessionList(slist{sli});
    SesList = {SesList(:).sessionName};

    for rli = 1:numel(rlist)        
        model = ['MTAC_BATCH-' fetSet ...
                 '_SR_' num2str(sampleRate) ...
                 '_NORM_' num2str(norm) ...         
                 '_REF_' rlist(rli).sessionName, '.' rlist(rli).mazeName '.' rlist(rli).trialName ...
                 '_STC_' rlist(rli).stcMode ...
                 '_NN_' num2str(nNeurons) ...
                 '_NI_' num2str(nIter) ...         
                 '_NN_multiPN_RAND_' rndMethod];
        load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'_PP',mapped,'.mat']));




        %figure
        hfig = figure(3923992);clf
        set(hfig,'Position',       [40,40,1500,500])
        set(hfig,'PaperPosition',  [0,0,1500/100,500/100])

        prop = 'accuracy';
        subplot(131);
        plot(cell2mat(cellfun(@subsref,ls, ...
                              repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
        xlim([0,5])
        ylim([0,100])
        ylabel(prop);
        title({['Training Set: ' rlist(rli).sessionName],...
               ['Randomization: ' rndMethod ],...
               ['Labeling Set: ' slist{sli}],...
               ['Feature  Set: ' fetSet]});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTick',1:numel(SesList));
        set(gca,'XTickLabel',SesList);
        set(gca,'XTickLabelRotation',90);
        pause(.1)

        prop = 'precision';
        subplot(132);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                                   repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'d-')
        xlim([0,7])
        hax = gca;
        hax.XTickLabelMode = 'manual';
        hax.XTickLabel = cat(2,{''},states,{''});
        ylim([0,100])
        ylabel(prop)
        title({['Training Set: ' rlist(rli).sessionName],...
               ['Randomization: ' rndMethod ],...
               ['Labeling Set: ' slist{sli}],...
               ['Feature  Set: ' fetSet]});
        legend(SesList,'location','southwest');
        pause(.1)

        prop = 'sensitivity';
        subplot(133);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                                   repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'d-')
        xlim([0,7])
        hax = gca;
        hax.XTickLabelMode = 'manual';
        hax.XTickLabel = cat(2,{''},states,{''});
        ylim([0,100])
        ylabel(prop)
        title({['Training Set: ' rlist(rli).sessionName],...
               ['Randomization: ' rndMethod ],...
               ['Labeling Set: ' slist{sli}],...
               ['Feature  Set: ' fetSet]});
        legend(SesList,'location','southwest');
        pause(.1)

        saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2/',...
                             ['Label_validation_WBSNT','_SL_',slist{sli}, ...
                            '_REF_',rlist(rli).sessionName,'.eps']),'epsc');

        
    end
end







figure,hold on,
fetSet  = 'fet_tsne_rev13';
k = 1;
htl = {};
set(gcf,'interpreter','none')
for sli = 1:numel(slist),
    for rli = 1:numel(rlist),    
        prop = 'accuracy';
        model = ['MTAC_BATCH-' fetSet ...
                 '_SR_' num2str(sampleRate) ...
                 '_NORM_' num2str(norm) ...         
                 '_REF_' rlist(rli).sessionName, '.' rlist(rli).mazeName '.' rlist(rli).trialName ...
                 '_STC_' rlist(rli).stcMode ...
                 '_NN_' num2str(nNeurons) ...
                 '_NI_' num2str(nIter) ...         
                 '_NN_multiPN_RAND_' rndMethod];

        ds = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'_PP',mapped,'.mat']));

        scat = scatter(k*ones([numel(ds.ls),1]),cell2mat(cellfun(@subsref,ds.ls, ...
                             repmat({substruct('.',prop)},[1,numel(ds.ls)]),'uniformoutput',false))'.*100,10);

        scat.MarkerFaceColor = scat.CData;
        k=k+1;
        htl = cat(2,htl,{['{' slist{sli} ' - ' rlist(rli).sessionName '}']});
    end
end

xlim([0.5,k-0.5]);
ylim([50,100])

set(gca,'XTickLabelMode','manual');
set(gca,'XTick',1:k-1);
set(gca,'XTickLabel',{});


set(gca,'XTickLabel',htl);
set(gca,'XTickLabelRotation',90);







