function fet = fet_respiration_reconstruction(Trial,varargin)
% [fet] = fet_respiration_reconstruction(Trial,varargin)
% 
% Uses feedforward neural network to reconstruct respiration signal from head-body movements.
% 
% Varargin:
%           NAME : TYPE    : DEFAULT VAL          : Description
%     sampleRate : Numeric : Trial.xyz.sampleRate : computational sampling rate
%        channel : Numeric : 2                    : ephys channel of respiratory signal
%            tag : String  : ''                   : save file tag 
%          train : Logical : false                : train feed forward neural network
%      overwrite : Logical : false                : flag for overwriting saved data
%

% DEFARGS ------------------------------------------------------------------------------------------
Trial = MTATrial.validate(Trial);
defargs = struct('sampleRate',   250,                                                            ...
                        'tag',   '',                                                             ...
                      'train',   false,                                                          ...
                  'overwrite',   false                                                           ...
);
[sampleRate,tag,train,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash(struct('function',mfilename()));
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,                                                   ... path
              [Trial.filebase,'.fet_respiration_reconstruction.',tag,'.mat'],... filename
              [],                                                            ... data
              sampleRate,                                                    ... sampleRate
              Trial.sync.copy,                                               ... sync object
              Trial.sync.data(1),                                            ... sync origin
              [],                                                            ... model
              'TimeSeries',                                                  ... data type
              [],                                                            ... extention
              'reconstructed respiration signal',                            ... name
              'fet_respiration_reconstruction',                              ... label
              'f');                                                            % key

if overwrite || ~exist(fet)
    
    pathNN = fullfile(Trial.path.project,'analysis',['nn_ff_respiration.mat']);
    
    twin = 64;
    
    pch = fet_HB_pitchB(Trial,sampleRate); 
    vxy = vel(filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',4,2.4,'low'),{'hcom'},[1,2]);
    ang = create(MTADang,Trial,filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',4,[30],'low'));
    rhm = fet_rhm(Trial,sampleRate);
    bresp = copy(ang);
    bresp.data = nunity(ang(:,'pelvis_root','spine_upper',3));
    bresp.filter('ButFilter',4,[4,15],'bandpass');
    lresp = copy(ang);
    lresp.data = nunity(ang(:,'pelvis_root','head_back',3));
    lresp.filter('ButFilter',4,[4,15],'bandpass');
    
    if train
        rfq = fet_respiration_freq(Trial,sampleRate,'overwrite',false);
        % DEFVAR lowpass filtered nasal pressure sensor
        ncpFilt = resample(filter(Trial.load('lfp',Trial.meta.channelGroup.respiration),...
                                  'ButFilter',4,16,'low'),rhm);
        ncpFilt.data = nunity(ncpFilt);

        gper = resample(cast([stc{'gper-groom'}],'TimeSeries'),rhm);
        gperSegs = gper.segs(1:twin/2:size(rhm,1),twin);
        gperSegs = all(gperSegs)';

        ncpSegs = ncpFilt.segs(1:twin/2:size(rhm,1),twin)';

        pchSegs = pch.segs(1:twin/2:size(rhm,1),twin); pchSegs = mean(pchSegs(:,:,1))';
        pchSegsMean = pch.segs(1:twin/2:size(rhm,1),twin); pchSegsMean = sq(mean(pchSegsMean(:,:,1)))';
        vxySegs = vxy.segs(1:twin/2:size(rhm,1),twin); vxySegs = sq(mean(vxySegs))';
        rfqSegs = rfq.segs(1:twin/2:size(rhm,1),twin); rfqSegs = sq(mean(rfqSegs))';

        rhmSegs   = clip(  rhm.segs(1:twin/2:size(rhm,1), twin)', -4, 4);
        lrespSegs = clip(lresp.segs(1:twin/2:size(rhm,1), twin)', -1, 1);
        brespSegs = clip(bresp.segs(1:twin/2:size(rhm,1), twin)', -1, 1);
        
        bind = nniz(rhmSegs)   & nniz(brespSegs) & nniz(lrespSegs) ...
               & nniz(ncpSegs) & nniz(pchSegs)   & nniz(vxySegs) & gperSegs;
        ncpSegs(~bind,:)     = [];
        brespSegs(~bind,:)   = [];
        lrespSegs(~bind,:)   = [];
        rhmSegs(~bind,:)     = [];
        pchSegsMean(~bind,:) = [];
        pchSegs(~bind,:)     = [];
        vxySegs(~bind,:)     = [];
        rfqSegs(~bind,:)     = [];

        clear('pchSegsInds','rfqSegsInds','vxySegsInds');
        vxySegsInds = discretize(vxySegs,logspace(-3,2,6));
        pchSegsInds(:,1) = discretize(pchSegsMean(:,1),linspace(-2,0.5,6));
        rfqSegsInds(:,1) = discretize(rfqSegs(:,1),linspace(0,15,6));

        dataTrain = [];
        dataTarget = [];

        for iter = 1:3
            sampleCount = 100;
            for rf = 1:5
                for v = 1:5
                    for ph = 1:5
                        ind = find(   v == vxySegsInds ...
                                      & ph == pchSegsInds(:,1) ...
                                      & rf == rfqSegsInds(:,1));
                        if numel(ind)>sampleCount
                            ind = ind(randsample(numel(ind),sampleCount,true));
                            dataTrain = cat(1,...
                                            dataTrain,...
                                            cat(2,...
                                                brespSegs(ind,:,:),...
                                                lrespSegs(ind,:,:),...
                                                rhmSegs(ind,:,:),...
                                                pchSegs(ind,:,:),...
                                                vxySegs(ind,:,:)) ...
                                            );
                            dataTarget = cat(1,dataTarget,ncpSegs(ind,:));
                        end
                    end
                end
            end
        end

        net = feedforwardnet([32,32]);
        net = train(net,dataTrain',dataTarget');
        save(pathNN,'net');
    else
        load(pathNN);
    end%if train
    
    alrespSegs = clip(lresp.segs(1:size(rhm,1),twin)',-1,1);
    abrespSegs = clip(bresp.segs(1:size(rhm,1),twin)',-1,1);
    arhmSegs = clip(rhm.segs(1:size(rhm,1),twin)',-4,4);
    apchSegs = pch.segs(1:size(rhm,1),twin);apchSegs = sq(mean(apchSegs(:,:,1)))';
    avxySegs = vxy.segs(1:size(rhm,1),twin);avxySegs = sq(mean(avxySegs))';
    
    nsegs = net(cat(2,abrespSegs,alrespSegs,arhmSegs,apchSegs,avxySegs)');
    nnsegs = (1./twin).*nsegs(1,:)';
    for s = 2:size(nsegs,1)
        nnsegs = nnsegs+(1./twin).*circshift(nsegs(s,:)',s);
    end    
    fet.data = nnsegs;
    
    % Use body feature for slow respiration
    sresp = copy(rhm);
    sresp.data = ang(:,'pelvis_root','spine_upper',3);
    sresp.filter('ButFilter',4,[1,6],'bandpass');

    ncpSynPow = copy(fet);
    ncpSynPow.data = sqrt(conv(nunity(fet.data).^2,gausswin(512)./sum(gausswin(512)),'same'));

    srespPeriods = ThreshCross(-ncpSynPow.data,-0.15,250);
    for s = 1:size(srespPeriods,1)
        fet.data(srespPeriods(s,1):srespPeriods(s,2)) =  ...
            sresp.data(srespPeriods(s,1):srespPeriods(s,2));
    end
    
    
    
% SAVE respiration frequency feature to file    
    fet.save();
else
% LOAD respiration feature from file    
    fet.load();
end


        