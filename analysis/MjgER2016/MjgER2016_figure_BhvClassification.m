% MjgER2016 Figure1
%
% Subplots:
%    A. Left, A schematic of the maze and the camera layout; 
%           center, image of the rat with tracking markers; 
%           Right, example of the marker skeleton with trajectories.
%    B. Top, An example feature matrix time series; 
%           middle, expert behavioral annotations; 
%           bottom, neural network behavioral annotations.
%    C. Labeling statistics of the neural network annotations with respect 
%           to the expert annotations. 
%           Left, overall accuracy; 
%           Middle, the precision of the state assignments; 
%           Right, the sensitivity of the state assignments. 
%    D. Joint probability distributions of two features which best separates one state 
%           from the others with countours of each state at the 95th percentile. 
%           Top left, sparation of the ?sit? behavior within the () and () space;
%           top right, separation of the ?rear? behavior within the () and () space; 
%           bottom left, separation of the ?groom? behavior within the () and () space; 
%           bottom right, ?walk?, ?turn?, and ?pause? within the () and () space.
%    E. Example traces relating the head pitch kinematics with respect to the respiration. 
%           From top to bottom, head pitch, bandpass-filtered (6-12Hz) head pitch, 
%           intra-nasal pressure, instantaneous respiration frequency, and ethogram.
%    F. Summary of the relationship of the head pitch dynamics to the respiration.
%           Coherence between the intra-nasal pressure and the head pitch angular velocity.
%           Diagram of the sniffing respiration cycle.
%           Head movements within the sniffing respiration cycle.
%    G. Top, the joint probability density; middle average intra-nasal pressure; 
%           bottom, average normalized head pitch, within the space of the phase difference 
%                   between the respiration and the head pitch angular velocity, and the 
%                   instantaneous head pitch frequency. 
%    H. Changes in the interaction between respiration and the head pitch angular velocity 
%           as a function of head pitch.
%    I. Schematic of the head during high and low pitch postures.
%    J. Probability distribution for normalized head pitch over multiple subjects; 
%           mean (black) ???confidence interval (grey).
%    K. Resultant length between pressure sensor and rhythmic head pitch related to normalized
%           head pitch, colored by phase difference.
%    L. Change in resultant length between high and low head posture in multiple subjects.
%    M. Change in pressure sensor frequency between high and low head posture in multiple subjects.
%    N. Probability distribution of head pich for each behavioral state.
%    O. Proportion of high and low head postures for each behavioral state.

% Supplementary Plots
% 
%    Behavioral state occpancy and transitions
%        For each session 
%            - accumulate the time spent in each state
%            - compute the transition matrix
%    Kinematics
%        
%
%    Time shifted Mutual Information of head and body speed demonstrats that the head movements lead
%        the body.
%
%    


global MTA_PROJECT_PATH
configure_default_args();


sessionListName = 'hand_labeled';
sessionList = get_session_list(sessionListName);
Trials = af(@(s) MTATrial.validate(s), sessionList);
numTrials = numel(Trials);

% STATE subspace data ----------------------------------------------------------
referenceTrial = 'jg05-20120317.cof.all';
states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_mis';
features = cf(@(T)  feval(featureSet,T,[],false),  Trials);
xyzs     = cf(@(T)  T.load('xyz'),                 Trials);
           cf(@(X,F)X.resample(F),                 xyzs,features);
Stc      = cf(@(T)  T.load('stc'),                 Trials);
stcm     = cf(@(S,X) stc2mat(S,X,states),          Stc, features);
cf(@(F,T) F.map_to_reference_session(T,referenceTrial),  features, Trials);
for s = 1:numTrials, features{s}.data(~nniz(xyzs{s}),:,:) = 0;end
nnizMat =  cat(1,cf(@(F) nniz(get(F,'data')), xyzs));
nnizMat = cat(1,nnizMat{:});
fetMat = cat(1,cf(@(F) get(F,'data'), features));
fetMat = cat(1,fetMat{:});
stcMat = cat(1,stcm{:});
%-------------------------------------------------------------------------------

sampleRate = 250;

% $$$ 
% $$$ Trial = Trials{3};
% $$$ xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD',sampleRate);
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ ncp = Trial.load('lfp',2);
% $$$ % $$$ ncp.resample(xyz);
% $$$ xts = [1:size(xyz,1)]./xyz.sampleRate;
% $$$ pch = copy(xyz);
% $$$ pch.data = ang(:,'hcom','nose',2);
% $$$ pch.filter('ButFilter',4,[4,15],'bandpass');
% $$$ rhm = fet_rhm(Trial,sampleRate);
% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ plot(xts,nunity(ang(:,'hcom','nose',2))*2+6.5)
% $$$ plot(xts,nunity(pch.data)+2,'b')
% $$$ plot(xts,-nunity(rhm.data)-1,'r')
% $$$ plot(xts,nunity(ncp.data)-4.5)
% $$$ xlim([860,865])


%featureMatrix = fet_mis(Trial);
% $$$ 
% $$$ states = {'rear','walk','turn','pause','groom','sit'};

%vxy = vel(filter(copy(preproc_xyz(Trials{tid},'SPLINE_SPINE_HEAD_EQD'),'ButFilter',4,2.4,'low'),{'bcom'},[1,2]);
% $$$ 
% $$$ % PLOT state jpdfs {normalizeSpineLength,normalizedSpineHight}, for each session
% $$$ figure(); 
% $$$ sv = {};
% $$$ sh = {};
% $$$ sax = reshape(tight_subplot(numStates,numTrials,[0.01,0.01],0.075,0.05),[numStates,numTrials])';
% $$$ for tid = 1:numel(Trials)
% $$$     stc = Trials{tid}.load('stc');
% $$$     [xyz,ss] = preproc_xyz(Trials{tid},'SPLINE_SPINE_HEAD_EQD');
% $$$     hvxy = fet_href_HXY(Trials{tid},[],[],'trb');
% $$$     sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
% $$$     sn = sum(sd(:,2:end-1),2)./sd(:,end);
% $$$     sv{tid} = copy(xyz);
% $$$     sv{tid}.data = sn;
% $$$     sd = mean(ss.data(:,:,3),2);
% $$$     sn = sd./mean(sd(hvxy(:,1)>10 & abs(hvxy(:,2))<20));
% $$$     sh{tid} = copy(xyz);
% $$$     sh{tid}.data = sn;
% $$$     for sts = 1:numStates
% $$$         axes(sax(sts,tid))
% $$$         hist2([log10(sv{tid}([stc{states{sts}}],1)),...
% $$$                log10(sh{tid}([stc{states{sts}}],1))],...
% $$$               linspace(-0.08,0.5,50),...
% $$$               linspace(-0.3,0.35,50));
% $$$         grid('on');
% $$$         set(gca(),'GridColor','w');
% $$$         caxis([0,prctile(nonzeros(get(get(gca,'Children'),'CData')),99)]);
% $$$         if sts~=numStates &&  tid~=numTrials
% $$$             sax(sts,tid).XTickLabel = {};
% $$$             sax(sts,tid).YTickLabel = {};            
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ % PLOT state jpdfs headVelocity {rostrocaudal,lateral}, for each session
% $$$ figure(); 
% $$$ sv = {};
% $$$ sh = {};
% $$$ sax = reshape(tight_subplot(numStates,numTrials,[0.01,0.01],0.075,0.05),[numStates,numTrials])';
% $$$ for tid = 1:numel(Trials)
% $$$     stc = Trials{tid}.load('stc');
% $$$     hvxy = fet_href_HXY(Trials{tid},[],[],'trb');
% $$$     for sts = 1:numStates
% $$$         axes(sax(sts,tid))
% $$$         hist2([hvxy([stc{states{sts}}],1),...
% $$$                hvxy([stc{states{sts}}],2)],...
% $$$               linspace(-20,80,50),...
% $$$               linspace(-60,60,50));
% $$$         grid('on');
% $$$         set(gca(),'GridColor','w');
% $$$         caxis([0,prctile(nonzeros(get(get(gca,'Children'),'CData')),99)]);
% $$$         if ~(sts==numStates &&  tid==1)
% $$$             sax(sts,tid).XTickLabel = {};
% $$$             sax(sts,tid).YTickLabel = {};            
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ % PLOT state jpdfs acom {height,velocity}, for each session
% $$$ figure(); 
% $$$ sax = reshape(tight_subplot(numStates,numTrials,[0.01,0.01],0.075,0.05),[numStates,numTrials])';
% $$$ for tid = 1:numel(Trials)
% $$$     stc = Trials{tid}.load('stc');
% $$$     xyz = filter(preproc_xyz(Trials{tid},'SPLINE_SPINE_HEAD_EQD'),'ButFilter',4,2.4,'low');
% $$$     ang = create(MTADang,Trials{tid},xyz);
% $$$     hvxy = fet_href_HXY(Trials{tid},[],[],'trb');    
% $$$     
% $$$     fhz = copy(xyz);
% $$$     fhz.data = xyz(:,'acom',3);
% $$$     fhz.data = fhz.data-mean(fhz(hvxy(:,1)>10 & abs(hvxy(:,2))<20));    
% $$$     fhdz = copy(xyz);
% $$$     fhdz.data = xyz(:,'acom',3);
% $$$     fhdz.data = (fhdz.data-circshift(fhdz.data,-1)).*fhdz.sampleRate./10;
% $$$     
% $$$     for sts = 1:numStates
% $$$         axes(sax(sts,tid))
% $$$         hist2([fhz([stc{states{sts}}],1),...
% $$$                fhdz([stc{states{sts}}],1)],...
% $$$               linspace(-40,100,50),...
% $$$               linspace(-20,20,50));
% $$$         grid('on');
% $$$         set(gca(),'GridColor','w');
% $$$         caxis([0,prctile(nonzeros(get(get(gca,'Children'),'CData')),99)]);
% $$$         if ~(sts==numStates &&  tid==1)
% $$$             sax(sts,tid).XTickLabel = {};
% $$$             sax(sts,tid).YTickLabel = {};            
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$     
% $$$     
% $$$ 
% $$$ 
% $$$ 
% $$$ fbxy = fet_bref_BXY(Trial,[],[],'SPLINE_SPINE_HEAD_EQD');
% $$$ 
% $$$ bvec = xyz(:,'hcom',[1,2])-xyz(:,'bcom',[1,2]);
% $$$ bvec = bsxfun(@rdivide,sq(bvec),sqrt(sum(bvec.^2,3)));
% $$$ avec = dot(sq(xyz(:,'acom',[1,2])-xyz(:,'bcom',[1,2])),bvec,2);
% $$$ avec = avec./mean(avec(vxy.data>10),'omitnan');
% $$$ 
% $$$ figure,
% $$$ subplot(411);
% $$$ hist2([ang(:,'bcom','hcom',2),avec],linspace(-0.5,pi/2,100),linspace(0,1.5,100));
% $$$ caxis([0,1000])
% $$$ subplot(412);
% $$$ ind = logical(get(cast([stc{'w'}],'TimeSeries'),'data'));
% $$$ hist2([ang(ind,'bcom','hcom',2),avec(ind)],linspace(-0.5,pi/2,100),linspace(0,1.5,100));
% $$$ caxis([0,1000])
% $$$ subplot(413);
% $$$ ind = logical(get(cast([stc{'p'}],'TimeSeries'),'data'));
% $$$ hist2([ang(ind,'bcom','hcom',2),avec(ind)],linspace(-0.5,pi/2,100),linspace(0,1.5,100));
% $$$ caxis([0,1000])
% $$$ subplot(414);
% $$$ ind = logical(get(cast([stc{'m'}],'TimeSeries'),'data'));
% $$$ hist2([ang(ind,'bcom','hcom',2),avec(ind)],linspace(-0.5,pi/2,100),linspace(0,1.5,100));
% $$$ caxis([0,1000])





% LOAD feature matrix example data
Trial = Trials{3}
stcHL = Trials{3}.load('stc');
stcML = Trials{3}.load('stc','msnn_ppsvd_raux');
fet = fet_mis(Trials{1});
xts = [1:size(fet,1)]./fet.sampleRate;


% LOAD labeling statistics
statsName = 'MTAC_STATS+TRN+hand_labeled+LBS+hand_labeled+fet_mis+SR10NN25NI100M1MREF+jg05-20120317.cof.all+N1NREF+hand_labeled+RNDWSBT+PRCT90+STS+wrnpms-multiSesPatNet';
%    stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));
stsOpt = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,'_pp.mat']));
states = stsOpt.states;
numStates = numel(states);
stateColors = [
               0.0000, 0.0000, 0.0000; ... walk
               0.3941, 0.0000, 0.6098; ... rear
               0.4235, 0.7333, 0.2352; ... turn    
               0.8139, 0.3861, 0.0000; ... pause
               0.8843, 0.2960, 0.2960; ... groom
               0.0000, 0.2607, 0.4843; ... sit
               ];

%startfig
[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],2,2,0.1,0.1);
globalXOffset = 0;
globalYOffset = 0;

% PLACEHOLDER for maze schematic and rat image and skelleton
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*4,                     ...
                              fig.subplot.height*1.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
text(sax(end),0.25,0.5,'IMAGE PLACEHOLDER');


% FEATURE matrix
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*4,                     ...
                              fig.subplot.height*2],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
imagesc([1:size(fet,1)]./fet.sampleRate,1:size(fet,2),nunity(fet.data)');
colormap('jet');
caxis([-3,3]);
xlim([430,460]);
sax(end).XTickLabel = {};

% EXPERT LABELS 
[yind, yOffSet, xind, xOffSet] = deal(5, 2.35, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*4,                     ...
                              fig.subplot.height*0.25],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
plotSTC(stcHL,1,'',states,stateColors,'staggeredStates',false);
xlim([430,460]);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
ylabel('EXP');

% EXPERT LABELS 
[yind, yOffSet, xind, xOffSet] = deal(6, 3.85, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*4,                     ...
                              fig.subplot.height*0.25],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
plotSTC(stcML,1,'',states,stateColors,'staggeredStates',false);
xlim([430,460]);
sax(end).YTickLabel = {};
ylabel('NN');
xlabel('Time (s)');

% STATE color legend
[yind, yOffSet, xind, xOffSet] = deal(4, 1, 5, -0.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*0.5,                    ...
                              fig.subplot.height*2],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
scatter(ones([numStates,1]),1:numStates,30,stateColors,'Filled')
sax(end).Visible= 'Off';
ylim([0,numStates+1])
title('Accuracy')
for sts = 1:numStates
text(1.3,sts-0.05,states{sts});
end
axis('ij');

% LABELING STATS
% ACCURACY 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 1, 5, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*0.5,                     ...
                              fig.subplot.height*1.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
    dat = cell2mat(af(@(x)x.accuracy,stsOpt.labelingStats))
plot(ones(size(dat)),dat.*100,'+')
ylim([83,87]);
title('Accuracy')


% PRECISION 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 1, 6, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');
dat = cell2mat(af(@(x)x.precision,stsOpt.labelingStats))';
xcoord = reshape((ones([numStates,1])*[1:numTrials])',[],1);
for sts = 1:numStates
    plot(xcoord(xcoord==sts),dat(xcoord==sts),'+','Color',stateColors(sts,:))
end
title('Precision')
xlim([0,numStates]+0.5);
sax(end).XTick = 1:numStates;
sax(end).XTickLabel = stsOpt.states;
sax(end).XTickLabelRotation = 90;
sax(end).YTick = 40:10:100;
grid(sax(end),'on');
ylim([40,100]);

% SENSITIVITY 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 1, 7, 1.25);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');

dat = cell2mat(af(@(x)x.sensitivity,stsOpt.labelingStats))';
xcoord = reshape((ones([numStates,1])*[1:numTrials])',[],1);
for sts = 1:numStates
    plot(xcoord(xcoord==sts),dat(xcoord==sts),'+','Color',stateColors(sts,:))
end
title('Sensitivity')
xlim([0,numStates]+0.5);
sax(end).YTickLabel = {};
sax(end).XTick = 1:numStates;
sax(end).XTickLabel = stsOpt.states;
sax(end).XTickLabelRotation = 90;
sax(end).YTick = 40:10:100;
sax(end).YTickLabel = {};
grid(sax(end),'on');
ylim([40,100]);


% LOAD DATA for respiration traces 
% ncp


configure_default_args();
meta = get_session_list_v2('ncp');
Trials = af(@(s) MTATrial.validate(s), meta);
numTrials = numel(Trials);


% Ed01-20140707.cof.all
% Ed01-20140709.cof.all
% Ed01-20140717.cof.all
% Ed05-20140528.cof.all
% Ed05-20140529.ont.all
% Ed10-20140815.cof.all
% Ed10-20140816.cof.all
% Ed10-20140817.cof.gnd


ncpStates = {'w+n+p',...
             'lloc+lpause','hloc+hpause',...
             'lloc','hloc','hpause','lpause','rear'};
nsts = numel(ncpStates);
titleText = {'fwdRHM (4-12) VS','Nasal Pressure Change'};
clear('xcomp','ycomp','zcomp','vdata');
tcomp = repmat(struct('data',[]),[1,nsts]);
xcomp = repmat(struct('data',[],'edgs',linspace([-pi,pi,15]),'label','Phase Difference','ctrs',[],'inds',[]),[1,nsts]);
ycomp = repmat(struct('data',[],'edgs',linspace([4,12,15]),'label',{'Nasal Presure Change','Frequenecy'},'ctrs',[],'inds',[]),[1,nsts]);
vdata = repmat(struct('data',[],'lim',[-1.2,0.1],'label','Head Pitch (rad)'),[1,nsts]);
zcomp = repmat(struct('count',[],'mean',[],'median',[],'std',[]),[1,nsts]);

% Manually select respiration troughs
% $$$ for tind = 1:numel(Trials),find_respiration_troughs(Trials{tind});end


sampleRate = 250;
for tind = 1:numel(Trials)
    Trial = Trials{tind};
    xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD',sampleRate);
    stcML = Trial.load('stc','msnn_ppsvd_raux');
    
    ncp = Trial.load('lfp',Trial.meta.channelGroup.respiration);
    ncpFilt = filter(copy(ncp),'ButFilter',4,[2,16],'bandpass');

    pch = fet_head_pitch(Trial,sampleRate);
    pch.map_to_reference_session(Trial,'Ed05-20140529.ont.all','');
    rhm = fet_rhm(Trial,sampleRate);
    rhmPhz = phase(rhm,[2,16]);

    %figure,plot(nunity(rhm.data));hold('on');plot(nunity(ncp.data));plot(nunity(rhm2.data));

    hvlf = fet_href_HXY(Trial,sampleRate,[],'trb');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    
    [amins,ncpSampleRate] = find_respiration_troughs(Trial);
    
    rmins = round((amins-1)/ncpSampleRate.*xyz.sampleRate);

% COMPUTE the interpolated phase based on the peaks detected in the ncp signal
    iphase = ncpFilt.copy();
    iphase.data(1:amins(1),1) = 0;
    for mm = 1:numel(amins)-1,
        iphase.data(amins(mm):amins(mm+1),1) = linspace(-pi,pi,amins(mm+1)-amins(mm)+1);
    end
    if amins(mm+1)<size(iphase,1)
        iphase.data(amins(mm+1):end,1) = 0;
    end
% COMPUTE the instantaneous respiratory frequency in rhm sampling rate
    ifreq = copy(pch);
    ifreq.data(1:rmins(1),1) = 0;
    for mm = 1:numel(rmins)-1,
        ifreq.data(rmins(mm):rmins(mm+1),1) = 1./((amins(mm+1)-amins(mm))./ncpSampleRate);
    end
    ifreq.data(rmins(mm+1):end,1) = 0;
    %ifreq = fet_respiration_freq(Trial,xyz.sampleRate,ncpChannel);
% COMPUTE the interpolated respiratory phase in rhm sampleRate
    riphase = copy(iphase);
    riphase.data = unwrap(iphase.data);
    riphase.resample(xyz);
    riphase.data = mod(riphase.data+2*pi,2*pi);
    riphase.data(riphase.data>pi) = riphase.data(riphase.data>pi)-2*pi;
    
    ncpRhmPDiff = circ_dist(riphase.data,rhmPhz.data);
    
% DOWNSAMPLE ncp to rhm
    rncp = resample(copy(ncp),xyz);
    
%% NCP x PCH -------------------------------------------------------------------
    for s = 1:nsts
    sind = logical(get(resample(cast([stcML{ncpStates{s}}],'TimeSeries'),xyz),'data'));
    srmins = rmins(sind(rmins));
    tcomp(s).data = cat(1,tcomp(s).data,tind.*ones([numel(srmins),1]));
    xcomp(s).data = cat(1,xcomp(s).data,ncpRhmPDiff(srmins,1));
    ycomp(s).data = cat(1,ycomp(s).data,ifreq(srmins,1));
    vdata(s).data = cat(1,vdata(s).data,pch(srmins,1));
    end
end


for s = 1:nsts
    xcomp(s).edgs = linspace([-pi,pi,32]);
    ycomp(s).edgs = linspace([4,12,22]);
    [xcomp(s),ycomp(s),zcomp(s)] = compute_2d_discrete_stats(xcomp(s),ycomp(s),vdata(s));
end


figure();
for s = 1:8
    subplot2(3,8,1,s); imagesc(xcomp(s).ctrs,ycomp(s).ctrs,zcomp(s).mean'); axis('xy');
        cax = colorbar(); ylabel(cax,['Mean ',vdata(s).label]); colormap('jet'); title(titleText); caxis([vdata(s).lim]);
    subplot2(3,8,2,s); imagesc(xcomp(s).ctrs,ycomp(s).ctrs,zcomp(s).std'); axis('xy');
        cax = colorbar(); ylabel(cax,['Std ',vdata(s).label]); ylabel(ycomp(s).label);
        colormap('jet');
    subplot2(3,8,3,s);
        imagesc(xcomp(s).ctrs,ycomp(s).ctrs,log10(zcomp(s).count)'); axis('xy');
        cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
        title(ncpStates{s});
        xlabel(xcomp(s).label);
        caxis([1,3.4]);
end

% jpdf
zflims = [0.0014,0.00325;...
          0.002,0.0044;...
          0.0015,0.00255;...
          0.0023,0.00515;...
          0.0016,0.00245;...
          0.0015,0.003125;...
          0.0016,0.00371;...
          0.0016,0.00285];

hcntr = {};
figure();

for s = 1:8
subplot2(3,8,1,s); imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,repmat(zcomp(s).mean,[2,1])'); axis('xy');
cax = colorbar(); ylabel(cax,['Mean ',vdata(s).label]); colormap('jet'); title(titleText); caxis([-1.25,1.25]);
subplot2(3,8,2,s); imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,repmat(zcomp(s).std,[2,1])'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',vdata(s).label]); ylabel(ycomp(s).label);
    colormap('jet');
subplot2(3,8,3,s);
    hold('on');
    zfilt = imgaussfilt(repmat(zcomp(s).count./sum(zcomp(s).count(:)),[2,1]),1.2);
    zfilt = zfilt./(0.5*sum(zfilt(:)));
    imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,zfilt'); axis('xy');
    xyp = cell([1,2]);
    [xyp{:}] = ndgrid([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs);
    [ccntr,hcntr{s}] = contour(xyp{:},zfilt,zflims(s,:).*[1,1]);   
    hcntr{s}.LineColor = 'k';hcntr{s}.LineWidth = 1;    
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet');     
    title(ncpStates{s});
    xlabel(xcomp(s).label);
    caxis([0,0.0085]);
end

figure
k = 1;
hcntr = {};
for s = [8,5,6,4,7]
    subplot(3,2,k);
    hold('on');
    zfilt = imgaussfilt(repmat(zcomp(s).count./sum(zcomp(s).count(:)),[2,1]),1.2);
    zfilt = zfilt./(0.5*sum(zfilt(:)));
    imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,zfilt'); axis('xy');
    xyp = cell([1,2]);
    [xyp{:}] = ndgrid([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs);
    [ccntr,hcntr{s}] = contour(xyp{:},zfilt,zflims(s,:).*[1,1]);   
    hcntr{s}.LineColor = 'k';hcntr{s}.LineWidth = 1;    
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet');     
    title(ncpStates{s});
    xlabel(xcomp(s).label);
    caxis([0,0.0085]);
    if k==1
        k=k+2;
    else
        k=k+1;
    end
end


figure
k = 1;
hcntr = {};
%for s = [8,5,4,6,7]
for s = [8,5,6,4,7]    
    subplot(5,1,k);
    hold('on');
    zfilt = imgaussfilt(repmat(zcomp(s).count./sum(zcomp(s).count(:)),[2,1]),1.2);
    zfilt = zfilt./(0.5*sum(zfilt(:)));
    imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,zfilt'); axis('xy');
    xyp = cell([1,2]);
    [xyp{:}] = ndgrid([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs);
    [ccntr,hcntr{s}] = contour(xyp{:},zfilt,zflims(s,:).*[1,1]);   
    hcntr{s}.LineColor = 'k';hcntr{s}.LineWidth = 1;    
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet');     
    title(ncpStates{s});
    xlabel(xcomp(s).label);
    caxis([0,0.0085]);
    k=k+1;
end



ccntr(:,ccntr(1,:)<2|ccntr(1,:)>6) = [];

inp = inpolygon(xyp{:},ccntr(1,:),ccntr(2,:));
%inp(15:26,:)=~inp(15:26,:);
%inp([19:21,],:)=~inp(19:26,:);
% $$$ inp(:,13:end) = false;
% $$$ inp([19:22],1:12) = true;
% $$$ inp([46:50],6:12) = true;
% $$$ inp([48:52],2:8) = true;
figure,imagesc(inp')
sum(zfilt(inp)./2)




% Spectral coherence and stuff

statesCoher = {'w+n+p+r&a','rear&a','hloc&a','hpause&a','lloc&a','lpause&a'};
coherS = {};
coherC = {};
sampleRate = 250;
pchBinEdgs = linspace([-pi/2,pi/2,21]);
pchBinCtrs = mean([pchBinEdgs(1:end-1);pchBinEdgs(2:end)]);

for tind = 1:numel(Trials)
    Trial = Trials{tind};
    xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD',sampleRate);
    ang =  create(MTADang,Trial,xyz);
    stcML = Trial.load('stc','msnn_ppsvd_raux');
% LOAD respiration signal (RESP)
    ncp = Trial.load('lfp',Trial.meta.channelGroup.respiration);
    resample(ncp,xyz);
% LOAD rhythmic head motion signal (RHM)
    rhm = fet_rhm(Trial,sampleRate);
    pch = fet_head_pitch(Trial,sampleRate);
    pch.map_to_reference_session(Trial,'Ed05-20140529.ont.all','');
    fpch = copy(pch);
    fpch.filter('ButFilter',4,[4,13],'bandpass');
% COMPUTE cross-spectrum of RESP and RHM
    specFet = copy(ncp);
    specFet.data = cat(2,specFet.data,rhm.data);
    specFet.data = cat(2,specFet.data,fpch.data);
    swOrd = 9;
    specArgs = struct('nFFT'     ,2^(swOrd),...
                   'Fs'       ,specFet.sampleRate,...
                   'WinLength',2^(swOrd-1),...
                   'nOverlap' ,2^(swOrd-1)*.5,...
                   'FreqRange',[.5,20]);
    [ys,fs,ts] = fet_spec(Trial,specFet,'mtcsdglong',false,[],specArgs,[],true);
% COMPUTE the coherence between RHM and RESP
    coherS{tind} = zeros([numel(fs),numel(statesCoher)]);
    coherP{tind} = zeros([numel(fs),10,numel(statesCoher)]);    
    for s = 1:numel(statesCoher)
        ind = [stcML{statesCoher{s}}];        
        for c = 2:3        
            coherS{tind}(:,s,c-1) = mean(abs(ys(ind,:,1,c)),'omitnan')./mean(sqrt(ys(ind,:,1,1).* ...
                                                              ys(ind,:,c,c)),'omitnan');
        end
        ypch = copy(pch);
        ypch.resample(ys);
        ypchInd = copy(ys);
        ypchInd.data = discretize(ypch.data,pchBinEdgs);
        ypchInd.data = ypchInd(ind);
        sys = ys(ind,:,:,:);
        for p = 1:numel(pchBinCtrs),
            pind = ypchInd.data==p;
            if sum(pind)>1
                coherP{tind}(:,p,s,c-1) = ...
                    mean(abs(sys(pind,:,1,c)),'omitnan')./mean(sqrt(sys(pind,:,1,1).*sys(pind,:,c,c)),'omitnan');
            end%if
        end%for p
    end%for s
end%for tind

coherN = mean(abs(ys(ind,:,1,3)),'omitnan')./mean(sqrt(ys(ind,:,1,1).*ys(ind,:,3,3)),'omitnan');

figure
k = 1;
hcntr = {};
%for s = [8,5,4,6,7]
statesCoher = {'w+n+p+r&a','rear&a','hloc&a','hpause&a','lloc&a','lpause&a'};
statesCoherLabels = {'All','Rear','HLoc','HPause','LLoc','LPause'};
sind = [2,3,4,5,6];
for s = [8,5,6,4,7]    
    subplot2(5,4,k,[1,2,3]);
    hold('on');
    zfilt = imgaussfilt(repmat(zcomp(s).count./sum(zcomp(s).count(:)),[2,1]),1.2);
    zfilt = zfilt./(0.5*sum(zfilt(:)));
    imagesc([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs,zfilt'); axis('xy');
    xyp = cell([1,2]);
    [xyp{:}] = ndgrid([xcomp(s).ctrs,xcomp(s).ctrs+2*pi],ycomp(s).ctrs);
    [ccntr,hcntr{s}] = contour(xyp{:},zfilt,zflims(s,:).*[1,1]);   
    hcntr{s}.LineColor = 'k';hcntr{s}.LineWidth = 1;    
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet');     
    title(ncpStates{s});
    caxis([0,0.0085]);
    xlim([pi/2,2*pi+pi/2])
    set(gca(),'XTick',[pi,3*pi/2,2*pi])
    set(gca(),'XTickLabel',{});
    if k==5,
        xlabel(xcomp(s).label);
        set(gca(),'XTickLabel',{'-\pi','-\pi/2','0'});
    end
    ylabel('Hz');
    subplot2(5,4,k,4);
    hold('on');
    cc = cf(@(c) c(:,sind(k)), coherS);
    cc = cat(2,cc{:});
    ts = tinv([0.025,0.975],numel(Trials)-1);
    CI =  bsxfun(@times,ts,sem');
    errorbar(mean(cc'),fs,CI(:,1)',CI(:,2)','horizontal')
    ylim([4,12]);
    xlim([0.4,0.75])
    if k==5,
        xlabel('Coherence')
    end
    k=k+1;    
end

% SUPFIG rhmXncp Coherence session x state ---------------------------------------------------------
figure();

clf();
for t = 1:8
    for s = 1:6
        subplot2(8,6,t,s);
        imagesc(pchBinCtrs,fs,coherP{t}(:,:,s));
        caxis([0.4,0.8]);
        if t ==1,
            title(statesCoherLabels{s});
        end
        if s==1,
            ylabel({Trials{t}.name,'Hz'});
        end
        colormap('jet');
        axis('xy');
    end
end
cax = colorbar()
ylabel(cax,'Coherence')
% end SUPFIG rhmXncp Coherence session x state -----------------------------------------------------




% SUPFIG rhmXncp spectrograms ----------------------------------------------------------------------
supfig = figure();
sax = tight_subplot(4,1,[0.03,0.03],[0.15,0.15]);        
tsNcp = [1:size(ncp,1)]./sampleRate;
cax = gobjects([4,1]);
supfig.CurrentAxes = sax(2);
    hold('on');
    imagesc(ts,fs,log10(ys(:,:,1,1)'));
    axis('xy');
    caxis([2,7.3]);
    colormap('jet');
    cax(2) = colorbar();
    axis('tight');    
    sax(2).XTickLabel = {};
    ylabel(supfig.CurrentAxes,{'RESP','Hz'});    
    ylabel(cax(2),{'Spectral','Power'});
supfig.CurrentAxes = sax(3);
    imagesc(ts,fs,log10(ys(:,:,2,2)'));
    axis('xy');
    caxis([-7,-3.5]);
    cax(3) = colorbar();
    axis('tight');    
    sax(3).XTickLabel = {}; 
    ylabel(supfig.CurrentAxes,{'RHM','Hz'});
    ylabel(cax(3),{'Spectral','Power'});
supfig.CurrentAxes = sax(4);    
    imagesc(ts,fs,(abs(ys(:,:,1,2))./(sqrt(ys(:,:,2,2)).*sqrt(ys(:,:,1,1))))');
    axis('xy');
    axis('tight');
    caxis([0,1]);
    cax(4) = colorbar();    
    linkaxes(sax(2:end),'x');
    xlim([10,60]);
    ylabel(supfig.CurrentAxes,{'RHM x RESP','Hz'});    
    xlabel(supfig.CurrentAxes,'Time (s)');    
    ylabel(cax(4),{'Coherence'});

supfig.CurrentAxes = sax(1);
    hold('on');
    plot(tsNcp,nunity(ncp(:,1)),'k');
    plot(tsNcp,nunity(rhm(:,1)).*2-4.75,'m');    
    xlim([33,37]);
    ylim([-8.5,4.5]);
    cax(1) = colorbar();
    supfig.CurrentAxes.YTickLabel = {};
    sax(1).XTickLabel = sax(1).XTick;    
    ylabel(supfig.CurrentAxes, 'RHM  RESP');
    title(supfig.CurrentAxes,{'Coherence between Rhythmic Head Motion and Respiration',''});
    
    grid(supfig.CurrentAxes,'on');
af(@(a) set(a,'Position',a.Position.*[0,1,0,1]+[0.15,0,0.7,0]), sax);
af(@(a) set(a,'Position',a.Position.*[0,1,0,1]+[0.87,0,0.02,0]), cax);
sax(1).Position(2) = sax(1).Position(2)+0.05;
sax(2).YTickLabel = sax(2).YTick;
fax = axes('Parent',supfig,'Position',[0,0,1,1],'Visible','off');
xlim(fax,[0,1]);
ylim(fax,[0,1]);

line(fax,...
     [sum(sax(2).Position([1,3]).*[1,0.46]),sax(2).Position(1)],...
     [sum(sax(2).Position([2,4])),sax(1).Position(2)-0.05],...
     'Color','k');
line(fax,...
     [sax(2).Position(1)].*[1,1],...
     [sax(1).Position(2)-[0.05,0.04]],...
     'Color','k');
line(fax,...
     [sum(sax(2).Position([1,3]))].*[1,1],...
     sax(1).Position(2)-[0.05,0.04],...
     'Color','k');

line(fax,...
     [sum(sax(2).Position([1,3]).*[1,0.54]),sum(sax(2).Position([1,3]))],...
     [sum(sax(2).Position([2,4])),sax(1).Position(2)-0.05],...
     'Color','k');
for p = 2:4,
patch(fax,...
     [sum(sax(p).Position([1,3]).*[1,0.46]),sum(sax(p).Position([1,3]).*[1,0.46]),...
      sum(sax(p).Position([1,3]).*[1,0.54]),sum(sax(p).Position([1,3]).*[1,0.54])],...
     [sax(p).Position([2]),sum(sax(p).Position([2,4])),...
      sum(sax(p).Position([2,4])),sax(p).Position([2])],...
      'k',...
         'FaceColor',0.95.*[1,1,1],...
         'FaceAlpha',0.2);
end
% end supfig rhmXncp spectrograms ------------------------------------------------------------------


    
figure,hold('on');
plot(nunity(rncp.data)),plot(ifreq.data)
    
figure();hold('on');
histogram(ifreq([stcML{'lloc+lpause'}]),linspace([1,15,30]),'Normalization','probability');
histogram(ifreq([stcML{'hloc+hpause'}]),linspace([1,15,30]),'Normalization','probability');

figure();hold('on');
histogram(ifreq([stcML{'lloc'}]),linspace([1,15,30]),'Normalization','probability');
histogram(ifreq([stcML{'hloc'}]),linspace([1,15,30]),'Normalization','probability');

figure();hold('on');
histogram(ifreq([stcML{'lpause'}]),linspace([1,15,30]),'Normalization','probability');
histogram(ifreq([stcML{'hpause'}]),linspace([1,15,30]),'Normalization','probability');

% $$$ acoher =mean(abs(ys(ind,:,1,2)),'omitnan')./mean(sqrt(ys(ind,:,1,1).*ys(ind,:,2,2)),'omitnan');
% $$$ % $$$ figure,plot(fs,acoher);
% $$$ hold('on'),plot(fs,acoher);

% coherence as a function of head pitch

ypch = copy(xyz);
ypch.data = ang(:,'hcom','nose',2);
ypch.resample(ys);
ypchInd = copy(ys);
ypchInd.data = discretize(ypch.data,linspace([-1.5,0.5,21]));
ypchInd.data = ypchInd(ind);

vxy = vel(filter(preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD'),'ButFilter',4,2.4,'low'),{'hcom'},[1,2]);
vxy.resample(ys);
ylvxy = copy(vxy);
ylvxy.data(ylvxy.data<1e-3) = 1e-3;
ylvxy.data = log10(ylvxy.data);
ylvxyInd = copy(ylvxy);
ylvxyInd.data = discretize(ylvxy.data,linspace([-2,2,21]));
ylvxyInd.data = ylvxyInd(ind);

sys = ys(ind,:,:,:);

pcoher = [];
for p = 1:20,
    pind = ypchInd.data==p;
    %pind = ypchInd.data==p & ylvxyInd.data>10;
    %pind = ylvxyInd.data==p;
    if sum(pind)>1
    pcoher(:,p) = mean(abs(sys(pind,:,1,2)),'omitnan')./mean(sqrt(sys(pind,:,1,1).*sys(pind,:,2,2)), ...
                                                      'omitnan');
    end
end
%figure,imagesc(linspace([-1.5,0.5,20]),fs,pcoher)
figure,imagesc(linspace([-2,2,20]),fs,pcoher)
axis('xy');
ylabel('Frequency');
xlabel('Head Pitch');
xlabel('log10 Head Speed');
title('hPch x NCP');
title('fRHM x NCP');





%% Feature for slow respiration 
rbm = copy(rhm);
rbm.data = nunity(minus(ang(:,'pelvis_root','spine_upper',3),ang(:,'pelvis_root','spine_middle',3)));
rbm.filter('ButFilter',4,[0.8,16],'bandpass');

rbm = copy(rhm);
rbm.data = minus(ang(:,'pelvis_root','spine_upper',3),ang(:,'pelvis_root','spine_middle',3));
rbm.data = rbm.data-medfilt1(rbm.data,251,'omitnan');
rbm.filter('ButFilter',4,[1,4],'bandpass');
rbm.data = nunity(rbm.data);
%rbm.data = nunity(rbm.data-medfilt1(rbm.data,251,'omitnan'));


rpm = copy(rhm);
rpm.data = minus(ang(:,'pelvis_root','spine_upper',3),ang(:,'pelvis_root','spine_middle',3));
rpm.data = rpm.data-medfilt1(rpm.data,251,'omitnan');
rpm.filter('ButFilter',4,[3,14],'bandpass');
rpm.data = nunity(rpm.data);
%rbm.data = nunity(rbm.data-medfilt1(rbm.data,251,'omitnan'));


rpm = copy(rhm);
rpm.data= nunity(circ_dist(ang(:,'pelvis_root','spine_upper',2),ang(:,'pelvis_root','spine_middle',2)));
rpm.filter('ButFilter',4,[0.8,16],'bandpass');

vxy = vel(filter(resample(copy(preproc_xyz(Trials{6},'SPLINE_SPINE_HEAD_EQD')),sampleRate),'ButFilter',4,2.4,'low'),{'bcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

figure,
subplot(211);
hold('on');
plot(xts,unity(rncp.data))
plot(xts,rbm.data*20)
plot(xts,vxy(:,1));
%plot(xts,rpm.data)
%plot(xts,nunity(rhm.data))
Lines([],-1.25,'k');
subplot(212);
plotSTC(stc,1,'',states,stateColors,'staggeredStates',false);
linkaxes(findobj(gcf(),'Type','Axes'),'x');



%/storage/eduardo/data/project/general/Ed01-20140707/Ed01-20140707.lfp
%/storage/eduardo/data/project/general/Ed01-20140709/Ed01-20140709.lfp
%/storage/eduardo/data/project/general/Ed01-20140717/Ed01-20140717.lfp
/storage/eduardo/data/project/general/Ed03-20140624/Ed03-20140624.lfp
/storage/eduardo/data/project/general/Ed03-20140625/Ed03-20140625.lfp
%/storage/eduardo/data/project/general/Ed05-20140529/Ed05-20140529.lfp
/storage/eduardo/data/project/general/Ed10-20140803/Ed10-20140803.lfp
/storage/eduardo/data/project/general/Ed10-20140804/Ed10-20140804.lfp
/storage/eduardo/data/project/general/Ed10-20140806/Ed10-20140807.lfp
/storage/eduardo/data/project/general/Ed10-20140808/Ed10-20140808.lfp
/storage/eduardo/data/project/general/Ed10-20140812/Ed10-20140812.lfp
/storage/eduardo/data/project/general/Ed10-20140814/Ed10-20140814.lfp
%/storage/eduardo/data/project/general/Ed10-20140815/Ed10-20140815.lfp
%/storage/eduardo/data/project/general/Ed10-20140816/Ed10-20140816.lfp
%/storage/eduardo/data/project/general/Ed10-20140817/Ed10-20140817.lfp
/storage/eduardo/data/project/general/Ed10-20140818/Ed10-20140818.lfp
/storage/eduardo/data/project/general/Ed10-20140820/Ed10-20140820.lfp
/storage/eduardo/data/project/general/Ed10-20140821/Ed10-20140821.lfp
/storage/eduardo/data/project/general/Ed10-20140822/Ed10-20140822.lfp
/storage/eduardo/data/project/general/Ed10-20140823/Ed10-20140823.lfp
/storage/eduardo/data/project/general/Ed10-20140825/Ed10-20140825.lfp
/storage/eduardo/data/project/general/Ed10-20140828/Ed10-20140828.lfp
/storage/eduardo/data/project/general/Ed10-20140829/Ed10-20140829.lfp
/storage/eduardo/data/project/general/Ed10-20140831/Ed10-20140831.lfp
/storage/eduardo/data/project/general/Ed10-20140901/Ed10-20140901.lfp
/storage/eduardo/data/project/general/Ed10-20140903/Ed10-20140903.lfp
/storage/eduardo/data/project/general/Ed10-20140905/Ed10-20140905.lfp
/storage/eduardo/data/project/general/Ed11-20150807/Ed11-20150807.lfp
/storage/eduardo/data/project/general/Ed11-20150811/Ed11-20150811.lfp
/storage/eduardo/data/project/general/Ed11-20150813/Ed11-20150813.lfp
/storage/eduardo/data/project/general/Ed11-20151006/Ed11-20151006.lfp
%%/storage/eduardo/data/project/general/Ed11-20151016/Ed11-20151016.lfp
/storage/eduardo/data/project/general/Ed11-20151020/Ed11-20151020.lfp
/storage/eduardo/data/project/general/Ed11-20151021/Ed11-20151021.lfp
/storage/eduardo/data/project/general/Ed11-20151104/Ed11-20151104.lfp
/storage/eduardo/data/project/general/Ed11-20151108/Ed11-20151108.lfp
/storage/eduardo/data/project/general/Ed12-20150807/Ed12-20150807.lfp
/storage/eduardo/data/project/general/Ed12-20150810/Ed12-20150810.lfp
/storage/eduardo/data/project/general/Ed12-20150811/Ed12-20150811.lfp
/storage/eduardo/data/project/general/Ed12-20150813/Ed12-20150813.lfp
%%/storage/eduardo/data/project/general/Ed12-20151014/Ed12-20151014.lfp
%%/storage/eduardo/data/project/general/Ed12-20151019/Ed12-20151019.lfp
%%/storage/eduardo/data/project/general/Ed12-20151020/Ed12-20151020.lfp
%%/storage/eduardo/data/project/general/Ed12-20151021/Ed12-20151021.lfp
%%/storage/eduardo/data/project/general/Ed12-20151104/Ed12-20151104.lfp

/storage/eduardo/data/processed/xyz/Ed12-20151104

ncpPchPDiff = circ_dist(ncpPhz.data,pchPhz.data);




ncpFreq = copy(ncpPhz);
ncpFreq.data = unwrap(ncpPhz.data);
%ncpFreq.filter('ButFilter',4,1,'low');
ncpFreq.data = 2.*pi./diff(ncpFreq.data);
figure,plot(ncpFreq.data)
ylim([0,20])

pchFreq = copy(pchPhz);
pchFreq.data = unwrap(pchPhz.data);
%pchFreq.filter('ButFilter',4,1,'low');
pchFreq.data = 2.*pi./diff(pchFreq.data);
figure,plot(pchFreq.data)
ylim([0,20])



[yind, yOffSet, xind, xOffSet] = deal(8, 0, 1,0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*8,                     ...
                              fig.subplot.height*3],                   ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');


shift = 10;
hold('on');
pt = [27.79,28.49;29.95,30.72;31.64,32.41];
for p = 1:size(pt,1)
   patch([pt(p,1),pt(p,1),pt(p,2),pt(p,2)],...
         [-10,20,20,-10]+shift,...
         [0.85,0.85,0.85],'EdgeColor','none');
end
plot(xts,nunity(ang(:,'hcom','nose',2))*3+9.5+shift,'k','LineWidth',1)
plot(xts,nunity(pch.data)*0.75+3+shift,'Color',[0.9,0.7,0],'LineWidth',1)
%plot(xts,nunity(rhm.data)-0.5,'r')
plot(lts,nunity(ncp.data).*0.75-1+shift,'b','LineWidth',1)
line(xts([1,end]),[1,1].*10./2-10+shift,'LineStyle','--','Color',[0.25,0.25,0.25])
line(xts([1,end]),[1,1].*4./2-10+shift,'LineStyle','--','Color',[0.25,0.25,0.25])
plot(xts,ifreq.data./2-10+shift,'Color',[0.3235, 0.6333, 0.1352],'LineWidth',1);
ylim([-10,14]+shift)
text(27.01,10./2-10+shift+0.5,'10Hz');
text(27.01,4./2-10+shift-0.5,'4Hz');
plotSTC(stc,1,'',states,stateColors,'staggeredStates',false);
xlim([27,33]);
sax(end).YTick = [0.5,3,4,9,10,12.5,13.5,20];
sax(end).YTickLabel = {'State','Frequency','Instantaneous','Sensor','Pressure','(5-12Hz)','Pitch','Pitch'};
xlabel('Time (s)');





% sit features --------------------------------------------------------------------
xfi = 14;  xlbl = 'Mean Body Speed log10(cm/s)';
xe = linspace(-3.5,2,20);
yfi = 9;   ylbl = 'Mid Body Height (mm)';
ye = linspace(50,115,20);
xc = mean([xe(1:end-1);xe(2:end)]);
yc = mean([ye(1:end-1);ye(2:end)]);
cords = cell([1,2]);
[cords{:}] = ndgrid(xc,yc);
stsC = {};
stsH = {};
thresh = {};
figure(9999);
for sts = 1:6
ind = nnizMat;
subplot2(3,6,1,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) ~= sts & nnizMat;
subplot2(3,6,2,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
subplot2(3,6,3,sts);hold('on');hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);
thresh{sts} = 1;
while sum(out(out>thresh{sts}))./sum(out(:))>0.90;
    thresh{sts} = thresh{sts}+1;
end
[stsC{sts},stsH{sts}] = contour(cords{:},out,[thresh{sts},thresh{sts}],'LineWidth',2,'Color',stateColors(sts,:));
end
figure(hfig);
[yind, yOffSet, xind, xOffSet] = deal(4, 1.5, 6, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);

hold('on');
ind = nnizMat;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);;
out(out<10) = nan;
set(pcolor(xe(1:end-1),ye(1:end-1),log10(out')),'EdgeColor','none');
for sts = 1:6
copyobj(stsH{sts},gca());
end
axis('tight');
caxis([0,4])
close(figure(9999));
% $$$ figure();
% $$$ for sts = 1:6;
% $$$ ind = stcMat(:,sts) == sts & nnizMat;
% $$$ clear('xcomp','ycomp','zcomp');
% $$$ vlabel = 'spine sinuosity (A.U.)';
% $$$ vdata = fetMat(ind,17); vClim = [1,2.5];
% $$$ titleText = '';
% $$$ xcomp.data = fetMat(ind,xfi); xcomp.label = xlbl;
% $$$ xcomp.edgs = xe;
% $$$ ycomp.data = fetMat(ind,yfi); ycomp.label = ylbl;
% $$$ ycomp.edgs = ye;
% $$$ [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vdata);
% $$$ zmean = zcomp.mean;   zmean(zcomp.count<thresh{sts}) = nan;
% $$$ zstd  = zcomp.std;    zstd(zcomp.count<thresh{sts}) = nan;
% $$$ zcount = zcomp.count; zcount(zcomp.count<thresh{sts}) = nan;
% $$$ subplot2(6,3,sts,1); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zmean'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet'); title(titleText); caxis([vClim]);
% $$$     if sts == 3, ylabel(ycomp.label);end
% $$$     grid('on');
% $$$ subplot2(6,3,sts,2); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zstd'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',vlabel]); caxis([0,0.9]);
% $$$     colormap('jet');
% $$$     grid('on');
% $$$ subplot2(6,3,sts,3);
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zcount'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     if sts == 6, xlabel(xcomp.label);end
% $$$     grid('on');    
% $$$ end




% rear features --------------------------------------------------------------------
%xfi = 14;  xlbl = 'Lower Body Speed log10(cm/s)';
xfi = 16;  xlbl = 'Upper Body Speed Z axis (cm/s)';
xe = linspace(-3,2,20);
yfi = 4;   ylbl = 'Upper Body Pitch (rad)';
ye = linspace(-1,pi/2,20);
xc = mean([xe(1:end-1);xe(2:end)]);
yc = mean([ye(1:end-1);ye(2:end)]);
cords = cell([1,2]);
[cords{:}] = ndgrid(xc,yc);
stsC = {};
stsH = {};
thresh = {};
nsts = stcMat(:,6) ~= 6;
figure(9999);
for sts = 1:5
ind = nnizMat;
subplot2(3,6,1,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) ~= sts & nnizMat & nsts;
subplot2(3,6,2,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
subplot2(3,6,3,sts);hold('on');hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);
thresh{sts} = 1;
while sum(out(out>thresh{sts}))./sum(out(:))>0.90;
    thresh{sts} = thresh{sts}+1;
end
[stsC{sts},stsH{sts}] = contour(cords{:},out,[thresh{sts},thresh{sts}],'LineWidth',2,'Color',stateColors(sts,:));
end
figure(hfig)
[yind, yOffSet, xind, xOffSet] = deal(4, 1.5, 7, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);

hold('on');
hold('on');
ind = nnizMat & nsts;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);;
out(out<10) = nan;
set(pcolor(xe(1:end-1),ye(1:end-1),log10(out')),'EdgeColor','none');
for sts = 1:5
copyobj(stsH{sts},gca());
end
axis('tight');
caxis([0,4]);
close(figure(9999));
% $$$ figure();
% $$$ for sts = 1:5;
% $$$ ind = stcMat(:,sts) == sts & nnizMat & nsts;
% $$$ clear('xcomp','ycomp','zcomp');
% $$$ vlabel = 'spine sinuosity (A.U.)';
% $$$ vdata = fetMat(ind,17); vClim = [1,2.5];
% $$$ titleText = '';
% $$$ xcomp.data = fetMat(ind,xfi); xcomp.label = xlbl;
% $$$ xcomp.edgs = xe;
% $$$ ycomp.data = fetMat(ind,yfi); ycomp.label = ylbl;
% $$$ ycomp.edgs = ye;
% $$$ [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vdata);
% $$$ zmean = zcomp.mean;   zmean(zcomp.count<thresh{sts}) = nan;
% $$$ zstd  = zcomp.std;    zstd(zcomp.count<thresh{sts}) = nan;
% $$$ zcount = zcomp.count; zcount(zcomp.count<thresh{sts}) = nan;
% $$$ subplot2(6,3,sts,1); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zmean'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet'); title(titleText); caxis([vClim]);
% $$$     if sts == 3, ylabel(ycomp.label);end
% $$$     grid('on');
% $$$ subplot2(6,3,sts,2); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zstd'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',vlabel]); caxis([0,0.9]);
% $$$     colormap('jet');
% $$$     grid('on');
% $$$ subplot2(6,3,sts,3);
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zcount'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     if sts == 6, xlabel(xcomp.label);end
% $$$     grid('on');    
% $$$ end



% groom features --------------------------------------------------------------------
xfi = 18;  xlbl = 'Spine Sinuosity (A.U.)';
xe = linspace(0,4,20);
yfi = 7;  ylbl = 'Lower Spine Height (mm)';
ye = linspace(0,60,20);
xc = mean([xe(1:end-1);xe(2:end)]);
yc = mean([ye(1:end-1);ye(2:end)]);
cords = cell([1,2]);
[cords{:}] = ndgrid(xc,yc);
stsC = {};
stsH = {};
thresh = {};
nsts = stcMat(:,6) ~= 6 & stcMat(:,2) ~= 2;
figure(9999);
for sts = [1,3:5]
ind = nnizMat;
subplot2(3,6,1,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) ~= sts & nnizMat & nsts;
subplot2(3,6,2,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
subplot2(3,6,3,sts);hold('on');hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);
thresh{sts} = 1;
while sum(out(out>thresh{sts}))./sum(out(:))>0.90;
    thresh{sts} = thresh{sts}+1;
end
[stsC{sts},stsH{sts}] = contour(cords{:},out,[thresh{sts},thresh{sts}],'LineWidth',2,'Color',stateColors(sts,:));
end
figure(hfig)
% INSERT SUBPLOT
[yind, yOffSet, xind, xOffSet] = deal(5, 0.5, 6, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
ind = nnizMat & nsts;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);;
out(out<10) = nan;
set(pcolor(xe(1:end-1),ye(1:end-1),log10(out')),'EdgeColor','none');
for sts = [1,3:5]
copyobj(stsH{sts},gca());
end
axis('tight');
caxis([0,4])
close(figure(9999));
% $$$ 
% $$$ figure();
% $$$ for sts = [1,3:5];
% $$$ ind = stcMat(:,sts) == sts & nnizMat & nsts;
% $$$ clear('xcomp','ycomp','zcomp');
% $$$ vlabel = 'spine sinuosity (A.U.)';
% $$$ vdata = fetMat(ind,17); vClim = [1,2.5];
% $$$ titleText = '';
% $$$ xcomp.data = fetMat(ind,xfi); xcomp.label = xlbl;
% $$$ xcomp.edgs = xe;
% $$$ ycomp.data = fetMat(ind,yfi); ycomp.label = ylbl;
% $$$ ycomp.edgs = ye;
% $$$ [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vdata);
% $$$ zmean = zcomp.mean;   zmean(zcomp.count<thresh{sts}) = nan;
% $$$ zstd  = zcomp.std;    zstd(zcomp.count<thresh{sts}) = nan;
% $$$ zcount = zcomp.count; zcount(zcomp.count<thresh{sts}) = nan;
% $$$ subplot2(6,3,sts,1); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zmean'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet'); title(titleText); caxis([vClim]);
% $$$     if sts == 3, ylabel(ycomp.label);end
% $$$     grid('on');
% $$$ subplot2(6,3,sts,2); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zstd'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',vlabel]); caxis([0,0.9]);
% $$$     colormap('jet');
% $$$     grid('on');
% $$$ subplot2(6,3,sts,3);
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zcount'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     if sts == 6, xlabel(xcomp.label);end
% $$$     grid('on');    
% $$$ end

% walk features --------------------------------------------------------------------
xfi = 6;  xlbl = 'Traj PPC (A.U.)';
xe = linspace(-0.3,1.1,20);
yfi = 11;  ylbl = 'Lower Spine Speed (cm/s)';
ye = linspace(-3,2,20);
xc = mean([xe(1:end-1);xe(2:end)]);
yc = mean([ye(1:end-1);ye(2:end)]);
cords = cell([1,2]);
[cords{:}] = ndgrid(xc,yc);
stsC = {};
stsH = {};
thresh = {};
nsts = stcMat(:,6) ~= 6 & stcMat(:,2) ~= 2 & stcMat(:,5) ~= 5;
figure(9999);
for sts = [1,3:4]
ind = nnizMat;
subplot2(3,6,1,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) ~= sts & nnizMat & nsts;
subplot2(3,6,2,sts);hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
subplot2(3,6,3,sts);hold('on');hist2(fetMat(ind,[xfi,yfi]),xe,ye);
ind = stcMat(:,sts) == sts & nnizMat;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);
thresh{sts} = 1;
while sum(out(out>thresh{sts}))./sum(out(:))>0.90;
    thresh{sts} = thresh{sts}+1;
end
[stsC{sts},stsH{sts}] = contour(cords{:},out,[thresh{sts},thresh{sts}],'LineWidth',2,'Color',stateColors(sts,:));
end
figure(hfig)
% INSERT SUBPLOT
[yind, yOffSet, xind, xOffSet] = deal(5, 0.5, 7, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*1,                     ...
                              fig.subplot.height*1],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
ind = nnizMat & nsts;
out = hist2(fetMat(ind,[xfi,yfi]),xe,ye);;
out(out<10) = nan;
set(pcolor(xe(1:end-1),ye(1:end-1),log10(out')),'EdgeColor','none');
for sts = [1,3:4]
copyobj(stsH{sts},gca());
end
axis('tight');
caxis([0,4])
colormap(jet);
close(figure(9999));

% $$$ figure();
% $$$ for sts = [1,3:4];
% $$$ ind = stcMat(:,sts) == sts & nnizMat & nsts;
% $$$ clear('xcomp','ycomp','zcomp');
% $$$ vlabel = 'spine sinuosity (A.U.)';
% $$$ vdata = fetMat(ind,17); vClim = [1,2.5];
% $$$ titleText = '';
% $$$ xcomp.data = fetMat(ind,xfi); xcomp.label = xlbl;
% $$$ xcomp.edgs = xe;
% $$$ ycomp.data = fetMat(ind,yfi); ycomp.label = ylbl;
% $$$ ycomp.edgs = ye;
% $$$ [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vdata);
% $$$ zmean = zcomp.mean;   zmean(zcomp.count<thresh{sts}) = nan;
% $$$ zstd  = zcomp.std;    zstd(zcomp.count<thresh{sts}) = nan;
% $$$ zcount = zcomp.count; zcount(zcomp.count<thresh{sts}) = nan;
% $$$ subplot2(6,3,sts,1); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zmean'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet'); title(titleText); caxis([vClim]);
% $$$     if sts == 3, ylabel(ycomp.label);end
% $$$     grid('on');
% $$$ subplot2(6,3,sts,2); 
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zstd'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',vlabel]); caxis([0,0.9]);
% $$$     colormap('jet');
% $$$     grid('on');
% $$$ subplot2(6,3,sts,3);
% $$$     set(pcolor(xcomp.edgs(1:end-1),ycomp.edgs(1:end-1),zcount'),'EdgeColor','none'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     if sts == 6, xlabel(xcomp.label);end
% $$$     grid('on');    
% $$$ end


rfrq = {};
bphz = {};
bvfs = {};
sampleRate = 250;
rstates = {'hloc','lloc'};
for tind = 1:numel(Trials)
    Trial = Trials{tind};
    xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD',sampleRate);
    stcML = Trial.load('stc','msnn_ppsvd_raux');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    pbl = phase(bvfl,[1,5]);
    [amins,ncpSampleRate] = find_respiration_troughs(Trial);
    rmins = round((amins-1)/ncpSampleRate.*xyz.sampleRate);
    % COMPUTE the instantaneous respiratory frequency in rhm sampling rate
    ifreq = copy(xyz);
    ifreq.data = zeros([size(ifreq,1),1]);
    ifreq.data(1:rmins(1),1) = 0;
    for mm = 1:numel(rmins)-1,
        ifreq.data(rmins(mm):rmins(mm+1),1) = 1./((amins(mm+1)-amins(mm))./ncpSampleRate);
    end
    ifreq.data(rmins(mm+1):end,1) = 0;
    for sts = 1:numel(rstates);
        hper = [stcML{rstates{sts},sampleRate}];
        hper.data(diff(hper.data,1,2)<200,:) = [];
        hmins = rmins(WithinRanges(rmins,hper.data));
        bvfs{tind,sts} = bvfl(hmins,1);
        bphz{tind,sts} = pbl(hmins,2);
        rfrq{tind,sts} = ifreq(hmins);
    end
end


hmins = hmins(bvfl(hmins,1)>30);
figure,hist2([pbl(hmins,2),ifreq(hmins)],linspace(-pi,pi,20),linspace(2,12,8));


figure,
ind = bvfs>40;
hist2([bphz(ind),rfrq(ind)],linspace(-pi,pi,8),linspace(2,12,8));

figure,
ind = rfrq>6 & rfrq<12;
hist2([bphz(ind),bvfs(ind)],linspace(-pi,pi,16),linspace(-10,80,8));


figure,
for tind = 1:numel(Trials)
    for sts = 1:numel(rstates)
    subplot2(numel(Trials),2,tind,sts);
    histcirc(bphz{tind,sts}(bvfs{tind,sts}>20&rfrq{tind,sts}>4),12);
    end
end

    

figure,histcirc(bphz(5000:6000),16)


