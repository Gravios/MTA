sli = 2
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


            xyz = Trial.load('xyz');
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
% $$$ ysm = MTADxyz('data',double(0<stc2mat(ds.stc{s},xyz,states)),'sampleRate',xyz.sampleRate); 
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


            %% COMPOSITE
            ang = create(MTADang,Trial,xyz);

            StcCor = ds.stc{s}.copy;
            rhh = [];
            rthresh = -0.2;
            tds = ds.d_state{s};
            [tpv,tps] = sort(ds.p_state{s},2,'descend');
            try
            for rp = ds.stc{s}{'r'}.data',
                rhh(end+1) = max(((xyz(rp',10,3)-refMean(4))./refStd(4)).*((ang(rp',3,4,2)-refMean(8))./refStd(8)));
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


            %% DISTANCE
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
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('w')}.data(wd<wthresh,:) = [];
            end



            %% Correct pause mislabels 
            %HEIGHT
            %StcCor = ds.stc{s}.copy;
            try
            pd = [];
            pthresh = -0.37;
            for rp = StcCor{'p'}.data',
                %for rp = StcHL{'p'}.data',
                pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = 6;%mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('p')}.data(pd<pthresh,:) = [];
            end



            %HEIGHT
            try
            pd = [];
            pthresh = -0.37;
            for rp = StcCor{'s'}.data',
                pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));    
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
            StcCor.states{StcCor.gsi('s')}.data(pd>pthresh,:) = [];
            end



            %% Correct sit Duration
            %% DURATION
            try
            pd = [];
            pthresh = 1.9;
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



            %% Correct turn Duration
            %% DURATION
            try
            pd = [];
            pthresh = 1.2;
            key = 'n';
            for rp = StcCor{key}.data',
                pd(end+1) = log10(abs(diff(rp)));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
            end

            
            %% DURATION
            try
            pppd = [];
            pthresh = 1.2;
            key = 'p';
            for rp = StcCor{key}.data',
                pd(end+1) = log10(abs(diff(rp)));
                if pd(end)<pthresh,
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
            end


            %% DURATION
            try
            pd = [];
            pthresh = 1.4;
            key = 'm';
            for rp = StcCor{key}.data',
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
            StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
            end


            %% HEIGHT
            %StcCor = ds.stc{s}.copy;
            try
            pd = [];
            pthresh = -0.37;
            for rp = StcCor{'p'}.data',
                %for rp = StcHL{'p'}.data',
                pd(end+1) = mean((xyz(rp',3,3)-refMean(2))./refStd(2));
                if pd(end)<pthresh
                    pind = rp(1):rp(2);
                    tds(pind,2) = 0;
                    tps(pind,:) = circshift(tps(pind,:),-1,2);
                    % maybe make a cat function for MTADepoch 
                    msts = 6;%mode(tps(pind,1));
                    
                    StcCor.states{StcCor.gsi(states{msts})}.data = ...
                        [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
                end
            end
            for sts = StcCor.states, sts{1}.clean; end
            StcCor.states{StcCor.gsi('p')}.data(pd<pthresh,:) = [];
            end
            


            %% DURATION
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
            errorEpochs = Trial.stc{'e'};
            if isempty(errorEpochs)
                errorEpochs = labelingEpochs;
            else
                errorEpochs.cast('TimeSeries');
            end
            
            
            %ind = cstc<0.05;
            ind = any(csm.data,2)&any(shl.data,2)&~errorEpochs.data;


            tcm = confmat(shl(ind&labelingEpochs.data,:),csm(ind&labelingEpochs.data,:)); % #DEP: netlab
            ls{s}.confusionMatrix = round(tcm./xyz.sampleRate,2);
            ls{s}.precision = round(diag(tcm)./sum(tcm,2),4).*100;
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







