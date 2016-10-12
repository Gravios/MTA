
%subjects = {'jg04','jg05','ER06'};
subjects = {'ER06_BHV','Ed10_BHV'};
subjects = {'jg05'};

stcMode = 'NN0317R';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
overwrite = false;


for subject = subjects
    subject = subject{1};

    sessionList = get_session_list(subject);
    for s = 1:numel(sessionList),    

        % label behavior
        Trial = MTATrial.validate(sessionList(s));

        try,  Trial.load('stc',stcMode); 
        catch err
            warning('MTAStateCollection not found, running label_bhv_nn.m')
            labelBhv_NN(Trial,...
                        stcMode,...
                        'jg05-20120317.cof.all',...
                        'hl_3_jg_r',...
                        [],[],[],[],[],[],[],{'loc','rear','pause','groom','sit'});
            Trial.load('stc',stcMode);
        end
        
        label_aux_bhv_reduced(Trial,stcMode,'overwrite',overwrite);        
        
        pft = pfs_2d_theta(Trial,[],[],overwrite);
        mrt = pft.maxRate;

        % Reduce clu list based on theta pfs max rate
        units = pft.data.clu(mrt>1);
         
        % Compute place fields and subsampled estimate
        for sts = 1:numel(states),
            defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
            defargs.units = units;
            defargs.states = states{sts};
            defargs = struct2varargin(defargs);        
            MTAAknnpfs_bs(Trial,defargs{:});      
        end
        
        % Parse place field features

        analysisFileName = fullfile(Trial.spath,[Trial.filebase,'_pfStats.mat']);        

        if ~exist(analysisFileName,'file') || overwrite,
            
            cluMap = pfkbs{1}.data.clu;

            for s = 1:numel(states),
                for u = 1:numel(units),
                    [pfkstats{s,u},pfkshuff{s,u},pfmstats{s,u}] = PlaceFieldStats(Trial,pfkbs{s},units(u));
                end
            end

            % Select the biggest baddest place field patch for all units
            for s = 1:numel(states),
                for k = 1:pfkbs{1}.parameters.numIter,
                    % Retrieve the patch center of mass from patch with the highest firing rate
                    pcom = ...
                    %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                    cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                            pfkshuff(s,:),...
                            repmat({k},[1,numel(units)]),...
                            'UniformOutput',false);

                    pind = ~cellfun(@isempty,pcom);
                    peakPatchCOM(s,k,pind,:) = ...
                        cell2mat(cellfun(@(x) x(:,1),...
                                         pcom(pind),...
                                         'uniformoutput',false))';

                    
                    % Retrieve the max patch Firing rate
                    peakPatchRate(s,k,:) = ...
                        sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                                   pfkshuff(s,:),...
                                   repmat({k},[1,numel(units)])));

                    
                    % Retrieve the patch area from patch with highest firing rate
                    parea = ...
                        cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                                pfkshuff(s,:),...
                                repmat({k},[1,numel(units)]),...
                                'UniformOutput',false);
                    pind = ~cellfun(@isempty,parea);
                    peakPatchArea(s,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                                  parea(pind),...
                                                                  'uniformoutput',false))');
                end
            end



            save(analysisFileName),'-v7.3',...
                'stcMode','states','cluMap',...
                'pfkshuff','pfkstats',...
                'peakPatchArea','peakPatchCOM','peakPatchRate');
        else
            load(analysisFileName);
        end

    end
end



%% Testing New Version of placefield stats
Trial = MTATrial.validate('jg05-20120310');
stcMode = 'NN0317R';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
Trial.load('stc',stcMode);

pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;
units = pft.data.clu(mrt>1);

for sts = 1:numel(states),
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end


cluMap = pfkbs{1}.data.clu;

figure,pfkbs{s}.plot(10);

for s = 1:numel(states),
    for u = 1:numel(units),
        [pfkstats{s,u},pfkshuff{s,u},pfmstats{s,u}] = PlaceFieldStats(Trial,pfkbs{s},units(u));
    end
end

for s = 1:numel(states),
    for k = 1:pfkbs{1}.parameters.numIter,
        % Retrieve the patch center of mass from patch with the highest firing rate
        pcom = ...
        %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
            cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                    pfkshuff(s,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);

        pind = ~cellfun(@isempty,pcom);
        peakPatchCOM(s,k,pind,:) = ...
            cell2mat(cellfun(@(x) x(:,1),...
                             pcom(pind),...
                             'uniformoutput',false))';

        
        % Retrieve the max patch Firing rate
        peakPatchRate(s,k,:) = ...
            sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                       pfkshuff(s,:),...
                       repmat({k},[1,numel(units)])));

        
        % Retrieve the patch area from patch with highest firing rate
        parea = ...
            cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                    pfkshuff(s,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);
        pind = ~cellfun(@isempty,parea);
        peakPatchArea(s,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                      parea(pind),...
                                                      'uniformoutput',false))');
    end
end



save(analysisFileName),'-v7.3',...
     'stcMode','states','cluMap',...
     'pfkshuff','pfkstats',...
     'peakPatchArea','peakPatchCOM','peakPatchRate');


load(analysisFileName);


Trial = MTATrial.validate('jg05-20120317');
Trial.load('stc','hl_3_jg_r');


xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',xyz);
[State, hmm, decode] = gausshmm(ang(ind.data,'head_back','head_front',2),2);


headPitchState = zeros([size(xyz,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;


figure,hold on
plot(ang(:,'head_back','head_front',2));
plot(headPitchState);


[decode] = hmmdecode(ang(:,'head_back','head_front',2),size(ang,1),hmm);


% Find most likely hidden state sequence using Viterbi method
State = decode(1).q_star;
figure,hold on
plot(ang(:,'head_back','head_front',2));
plot(decode(1).q_star);



figure,hold on
plot(ang(:,'head_back','head_front',2));
plot(headPitchState);


figure,hold on
eds = linspace(-pi/2,pi/2,200);
ind = headPitchState(:,1)==1;
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
ind = headPitchState(:,2)==0;
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'r';
hs.EdgeColor = 'r';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;


% HMM head pitch during locomotion and pause states
angThreshHMM = mean([prctile(ang(headPitchState(:,1)==1,'head_back','head_front',2),2),prctile(ang(headPitchState(:,2)==1,'head_back','head_front',2),98)]);


% MEAN rhm ang threshold
compute_session_rhm_distribution(Trial,[],'loc+pause','hl_3_jg_r');

afig = hgload(fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig'));
rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','compute_session_rhm_distribution-hangle'));
rhm_distrb = rhm_distrb(1);
figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
mrhmp(mrhmp==0) = nan;
rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
    +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));




s = 1;
Trial = MTATrial.validate(sessionList(s));
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
ang = create(MTADang,Trial,xyz);


figure,hold on
eds = linspace(-pi/2,pi/2,200);
ind = [Trial.stc{'lloc'}];
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
ind = [Trial.stc{'hloc'}];
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'r';
hs.EdgeColor = 'r';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;






[rhm,fs] = fet_rhm(Trial,[],'mtchglong');
trhm = rhm.copy;
rhm.data  = log10(rhm.data);
rhm.data(rhm<-9) = nan;
rhm.data(nniz(rhm.data))=nan;
vel = xyz.vel(1,[1,2]);
vel.resample(rhm);
vnn = nniz(vel);
rhm.data = (rhm.data-repmat(nanmean(rhm(vnn,:)),[rhm.size(1),1]))...
           ./repmat(nanstd(rhm(vnn,:)),[rhm.size(1),1]);

rhmp = rhm.copy;
rhmp.data = nanmean(rhm(:,6<fs&fs<12),2);

ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',rhmp);
ind.data = ind.data&nniz(rhmp);
[State, hmm, decode] = gausshmm(rhmp(ind.data),2);



headPitchState = zeros([size(rhmp,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;
headPitchState(ind.data==1,2) = State'~=stsInd;


figure,hold on
plot(rhmp.data);
plot(headPitchState(:,1));
plot(headPitchState(:,2));

figure,hold on
plot(ang(headPitchState(:,1)==1,'head_back','head_front',2));
plot(headPitchState);

hps = rhmp.copy;
hps.data = headPitchState+eps;
hps.resample(xyz);

figure,hold on
eds = linspace(-pi/2,pi/2,200);
sind = hps(:,1)>0.5;
hs = bar(eds,histc(ang(sind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
sind = hps(:,2)>0.5;
hs = bar(eds,histc(ang(sind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'r';
hs.EdgeColor = 'r';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;


% HMM head pitch during locomotion and pause states
rhmThreshHMM = mean([prctile(ang(hps(:,1)>0.5,'head_back','head_front',2),98),prctile(ang(hps(:,2)>0.5,'head_back','head_front',2),2)]);








