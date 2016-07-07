function stc = optimize_stc_transition_mis_single_session(Trial,Stc,varargin)

% DEFARGS ----------------------------------------------------------------------
defargs = struct('RefTrial',  'jg05-20120317.cof.all',                       ...
                 'refStcMode','hand_labeled_rev3_jg',                        ...
                 'fetSet', 'fet_mis',                                        ...
                 'stcMode', '',                                              ...
                 'tag_preprocessing', '+seh+',                               ...
                 'tag_postprocessing','_PP',                                 ...
                 'sampleRate', 12,                                           ...
                 'nNeurons',   100,                                          ...
                 'nIter',      100,                                          ...
                 'states',    {{'walk','rear','turn','pause','groom','sit'}},...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',       true,                                         ...
                 'map2ref',       true                                          ...
);
%-----------------------------------------------------------------------------

[RefTrial,refStcMode,fetSet,stcMode,tag_preprocessing,tag_postprocessing,sampleRate,...
 nNeurons,nIter,states,rndMethod,norm,map2ref] = DefaultArgs(varargin,defargs,'--struct');

RefTrial = MTATrial.validate(RefTrial);

% MAIN -------------------------------------------------------------------------

model = ['MTAC_BATCH-' tag_preprocessing fetSet ...
         '_SR_' num2str(sampleRate) ...
         '_NORM_' num2str(norm) ...         
         '_REF_' RefTrial.filebase ...
         '_STC_' refStcMode ...
         '_NN_' num2str(nNeurons) ...
         '_NI_' num2str(nIter) ...         
         '_NN_multiPN_RAND_' rndMethod];

RefTrial.load('stc',refStcMode);

rfet = feval(fetSet,RefTrial,sampleRate,false);
[~,refMean,refStd] = nunity(rfet(RefTrial.stc{'a'},:));


Trial = MTATrial.validate(Trial);


if map2ref,
    mapped = '-map2ref';
else
    mapped = '';
end

ds = load(fullfile(Trial.spath,[Trial.filebase,'.','labelBhv_NN','-',model,mapped,'.mat']));

%% Pre-process features

features = feval(fetSet,Trial,sampleRate,false);
features.map_to_reference_session(Trial,RefTrial);
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


xyz = Trial.load('xyz','seh');
ss = Trial.load('fet','3dssh');
ss.resample(xyz);


xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);

ang = create(MTADang,Trial,xyz);

if ~isempty(stcMode),
    StcHL = Trial.load('stc',stcMode);
    shl = MTADxyz('data',double(~~stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
end

%% COMPOSITE: rear
StcCor = ds.stc;
rhh = [];
rthresh = 20;
tds = ds.d_state;
[tpv,tps] = sort(ds.p_state,2,'descend');
try
    for rp = StcCor{'r'}.data',
        %for rp = StcHL{'r'}.data',
        rhh(end+1) = max(xyz(rp',10,3).*ang(rp',3,4,2));
        %rhh(end+1) = min(((xyz(rp',10,3)-refMean(4))./refStd(4)).*((ang(rp',3,4,2)-refMean(9))./refStd(9)));
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


%% Correct Groom Duration
%% DURATION
try
    pd = [];
    pthresh = 2;
    for rp = StcCor{'m'}.data',
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
    StcCor.states{StcCor.gsi('m')}.data(pd<pthresh,:) = [];
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
    pthresh = 1.6;
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



ls = [];
if ~isempty(stcMode),
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
    ls.confusionMatrix = round(tcm./xyz.sampleRate,2);
    ls.precision = round(diag(tcm)./sum(tcm,2),4)'.*100;
    ls.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
    ls.accuracy = sum(diag(tcm))/sum(tcm(:));

end

StcCor.updateMode([StcCor.mode,tag_postprocessing]);
stc = StcCor.copy;
stc.save(1);


mapped = '-map2ref';

save(fullfile(Trial.spath,...
              [Trial.filebase,'.',mfilename,'-',model,tag_postprocessing,mapped,'.mat']),...
     '-v7.3','nNeurons','nIter','sampleRate','model',...
     'fetSet','rndMethod','states','stc','ls');





