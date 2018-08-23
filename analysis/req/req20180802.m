% phase precession as function of {drz, body pitch, head pitch}  

MjgER2016_load_data();

MjgER2016_load_bhv_erpPCA_scores();

[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
    req20180123_ver5(Trials,[],[],false,false)


trialIndex = 20;


% LOAD Trial
% LOAD state collection
% GET unit subset
% LOAD position (xyz) data
% COMPUTE intermarker angles and distances
Trial = Trials{trialIndex};
stc = Trial.load('stc','msnn_ppsvd_raux');
unitSubset = units{trialIndex};
xyz = preproc_xyz(Trial,'trb');  if xyz.sampleRate > 120,  xyz.resample(120);  end
ang = create(MTADang,Trial,copy(xyz));
fet = fet_HB_pitchB(Trial);
%hrf = fet_href_H(Trial);

states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};


% LOAD theta behavioral fields
% LOAD theta place fields
% COMPUTE direction rate zone of each unit based on theta place fields
pfb = pfd{trialIndex,1};
pft = pfs_2d_theta(Trial,unitSubset);
pfs = pfs_2d_states(Trial,unitSubset,[],states);
drz = MTADfet.encapsulate(Trial,                                                 ... MTATrial object
                          compute_drz(Trial,unitSubset,pft,'feature',xyz.copy()),... data 
                          xyz.sampleRate,                                        ... sample rate
                          'directional rate zones',                              ... name
                          'drz',                                                 ... label
                          'z');%                                                     key
ddz = MTADfet.encapsulate(Trial,                                                 ... MTATrial object
                          compute_ddz(Trial,unitSubset,pft,'feature',xyz.copy()),... data 
                          xyz.sampleRate,                                        ... sample rate
                          'directional distance zones',                          ... name
                          'ddz',                                                 ... label
                          'd');%                                                     key
% LOAD local field potential (LFP)
% COMPUTE phase between 5-12 Hz from LFP
% LOAD spike time objects                          
%lfp = Trial.load('lfp',[68,72,76,82]);
lfp = Trial.load('lfp',[8]);                          
phz = phase(resample(copy(lfp),xyz),[5,12]);
spk = Trial.load('spk',xyz.sampleRate,[],unitSubset,'deburst');



statePlotParam = {};

% THETA 
statePlotParam{end+1}.state = stc{'theta-groom-sit',xyz.sampleRate}; 
statePlotParam{end}.label = 'theta';
statePlotParam{end}.bpLims = [-0.7,1.3];
statePlotParam{end}.hpLims = [-1.1,1.3];
statePlotParam{end}.hbpLims = [-1.1,1.3];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,250];
% REAR 
statePlotParam{end+1}.state = stc{'rear&theta',xyz.sampleRate}; 
statePlotParam{end}.label = 'rear';
statePlotParam{end}.bpLims = [-0.6,1.3];
statePlotParam{end}.hpLims = [-1.1,1.3];
statePlotParam{end}.hbpLims = [-1.1,1.3];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,250];
% HLOC 
statePlotParam{end+1}.state = stc{'hloc&theta',xyz.sampleRate}; 
statePlotParam{end}.label = 'high loc';
statePlotParam{end}.bpLims = [-0.7,0.1];
statePlotParam{end}.hpLims = [-0.5,0.7];
statePlotParam{end}.hbpLims = [-0.5,0.7];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,100];
% HPAUSE 
statePlotParam{end+1}.state = stc{'hpause&theta',xyz.sampleRate}; 
statePlotParam{end}.label = 'high pause';
statePlotParam{end}.bpLims = [-0.7,0.1];
statePlotParam{end}.hpLims = [-0.5,0.7];
statePlotParam{end}.hbpLims = [-0.5,0.7];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,110];
% LLOC 
statePlotParam{end+1}.state = stc{'lloc&theta',xyz.sampleRate}; 
statePlotParam{end}.label = 'low loc';
statePlotParam{end}.bpLims = [-0.7,0];
statePlotParam{end}.hpLims = [-1.1,-0.1];
statePlotParam{end}.hbpLims = [-1.1,0];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,70];
% LPAUSE 
statePlotParam{end+1}.state = stc{'lpause&theta',xyz.sampleRate}; 
statePlotParam{end}.label = 'low pause';
statePlotParam{end}.bpLims = [-0.7,-0.1];
statePlotParam{end}.hpLims = [-1.1,-0.1];
statePlotParam{end}.hbpLims = [-1.1,0];
statePlotParam{end}.hbdLims = [-pi/2,pi/2];
statePlotParam{end}.zLims = [20,70];

% PARAMATIZATION of statePlotParams
% $$$ x = 1;
% $$$ munits = unitSubset(x:x+9);
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position = [1,1,60,20];
% $$$ s = 1;
% $$$ for u = 1:numel(munits),
% $$$     res = spk(munits(u));
% $$$     res = res(WithinRanges(res,statePlotParam{s}.state.data));
% $$$     subplot2(4,10,1,u);
% $$$     plot(repmat(ang(res,'spine_middle','spine_upper',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     xlim(statePlotParam{s}.bpLims);
% $$$     ylim([-pi,pi*3]);
% $$$     subplot2(4,10,2,u);
% $$$     plot(repmat(ang(res,'head_back','head_front',2),[2,1]), [phz(res,1);phz(res,1)+2*pi],'.');
% $$$     xlim(statePlotParam{s}.hpLims);
% $$$     ylim([-pi,pi*3]);
% $$$     subplot2(4,10,3,u);    
% $$$     plot(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')
% $$$     xlim([-1,1])
% $$$     ylim([-pi,pi*3]);
% $$$     subplot2(4,10,4,u);
% $$$     plot(repmat(xyz(res,'hcom',3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     xlim(statePlotParam{s}.zLims);
% $$$     ylim([-pi,pi*3]);
% $$$ end





%function [P,phzStats,Rmax] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units)
%
% Compute phase precession coefficients between distance restricted drz and lfp phase
% 
% Varargin:
%
%    drz:        directional rate zones
%    ddz:        directional distance zones
%    phz:        Local field potential phase
%    spk:        MTASpk object which holds spike events
%    units:      list of units for computation
%    distThresh: limit phase precession analysis to specified radius around placefield center
%    mResults:   number of best fits to be returned

dthresh = 250;
clear('ppHH','ppBP','ppHB','ppHP','ppDRZ');
for unit = 1:numel(unitSubset),
    for s = 1:numel(statePlotParam),
        res = spk(unitSubset(unit));
        res = res(WithinRanges(res,statePlotParam{s}.state.data));
        if ~isempty(res),
            [ppHH.P(unit,:,s),ppHH.phzStats(unit,:,s),ppHH.Rmax(unit,s)] = ...
                MjgER2016_phasePrecession(Trial,xyz(res,'hcom',3),...
                                          ddz(res,unit),phz(res,1),[],unit,dthresh,1);
            [ppBP.P(unit,:,s),ppBP.phzStats(unit,:,s),ppBP.Rmax(unit,:,s)] = ...
                MjgER2016_phasePrecession(Trial,fet(res,2),...
                                          ddz(res,unit),phz(res,1),[],unit,dthresh,1);
            [ppHB.P(unit,:,s),ppHB.phzStats(unit,:,s),ppHB.Rmax(unit,:,s)] = ...
                MjgER2016_phasePrecession(Trial,fet(res,1),...
                                          ddz(res,unit),phz(res,1),[],unit,dthresh,1);
            [ppHP.P(unit,:,s),ppHP.phzStats(unit,:,s),ppHP.Rmax(unit,:,s)] = ...
                MjgER2016_phasePrecession(Trial,ang(res,'head_back','head_front',2),...
                                          ddz(res,unit),phz(res,1),[],unit,dthresh,1);
            [ppDRZ.P(unit,:,s),ppDRZ.phzStats(unit,:,s),ppDRZ.Rmax(unit,:,s)] = ...
                MjgER2016_phasePrecession(Trial,drz(res,unit),...
                                          ddz(res,unit),phz(res,1),[],unit,dthresh,1);
        end
    end
end





hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,60,30];
ny = 6;
nx = 10;


for unit = unitSubset,
clf(hfig);
pfsMaxRate = max(cell2mat(cf(@(p) p.maxRate(unit),pfs)));

subplot2(ny,nx,2,1);
plot(pfb,unit,'mean',true,[0,pfsMaxRate],false,'colorMap',@jet);
cax = colorbar(gca);
cax.Position(1) = cax.Position(1)+0.04;    

uind = find(unit==unitSubset);

for s = 1:numel(statePlotParam),
    res = spk(unit);
    res = res(WithinRanges(res,statePlotParam{s}.state.data));
    res(abs(ddz(res,unit==unitSubset))>250) = [];
    res(drz(res,unit==unitSubset)>0) = [];    
    
    if isempty(res), continue, end

    p = 2;    
    % Place field
    subplot2(ny,nx,s,p);p=p+2;
    plot(pfs{s},unit,'mean',true,[0,pfsMaxRate],true,'colorMap',@jet);
    cax = colorbar(gca);
    cax.Position(1) = cax.Position(1)+0.04;
    title(statePlotParam{s}.label);
    
% $$$     % height    
% $$$     subplot2(ny,nx,s,p);p=p+1;
% $$$     plot(repmat(xyz(res,'hcom',3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     xlim(statePlotParam{s}.zLims);
% $$$     ylim([-pi,pi*3]);
% $$$     ylabel('Theta Phase');        
% $$$     if s==1,                     title('Height'); end
% $$$     if s==numel(statePlotParam), xlabel('mm');    end
% $$$ 
% $$$     % body pitch    
% $$$     subplot2(ny,nx,s,p);p=p+1;
% $$$     plot(repmat(ang(res,'spine_middle','spine_upper',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     xlim(statePlotParam{s}.bpLims);
% $$$     ylim([-pi,pi*3]);
% $$$     if s==1,                  title('Body Pitch');    end    
% $$$     if s==numel(statePlotParam), xlabel('rad');       end    
% $$$     
% $$$     % head pitch
% $$$     subplot2(ny,nx,s,p);p=p+1;
% $$$     plot(repmat(ang(res,'head_back','head_front',2),[2,1]), [phz(res,1);phz(res,1)+2*pi],'.');
% $$$     xlim(statePlotParam{s}.hpLims);
% $$$     ylim([-pi,pi*3]);
% $$$     if s==1,                     title('Head Pitch'); end        
% $$$     if s==numel(statePlotParam), xlabel('rad');       end        
% $$$     
% $$$     % drz
% $$$     subplot2(ny,nx,s,p);p=p+1;
% $$$     plot(repmat(drz(res,unit==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')
% $$$     xlim([-1,1])
% $$$     ylim([-pi,pi*3]);
% $$$     if s==1,                     title('DRZ');      end   
% $$$     if s==numel(statePlotParam), xlabel('AU');  end        
% $$$ 
% $$$ 
% $$$     
    % height    
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(xyz(res,'hcom',3),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace([statePlotParam{s}.zLims,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1),5,1)');    
    axis('xy');    
    ylabel('Theta Phase');        
    if s==1,                     title('Height'); end
    if s==numel(statePlotParam), xlabel('mm');    end
    colormap(gca,'default');    
    hold('on');
    plot(statePlotParam{s}.zLims,2*pi*ppHH.P(uind,1,s)*statePlotParam{s}.zLims+ppHH.P(uind,2,s),'-m','LineWidth',1);
    plot(statePlotParam{s}.zLims,2*pi*ppHH.P(uind,1,s)*statePlotParam{s}.zLims+ppHH.P(uind,2,s)+2*pi,'-m','LineWidth',1)
    

    % body pitch    
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(fet(res,2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace([statePlotParam{s}.bpLims,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1)',5,1));    
    axis('xy');    
    if s==1,                  title('Body Pitch');    end    
    if s==numel(statePlotParam), xlabel('rad');       end    
    colormap(gca,'default');    
    hold('on');
    plot(statePlotParam{s}.bpLims,2*pi*ppBP.P(uind,1,s)*statePlotParam{s}.bpLims+ppBP.P(uind,2,s),'-m','LineWidth',1);
    plot(statePlotParam{s}.bpLims,2*pi*ppBP.P(uind,1,s)*statePlotParam{s}.bpLims+ppBP.P(uind,2,s)+2*pi,'-m','LineWidth',1)
    
    
    % head pitch
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(ang(res,'head_back','head_front',2),[2,1]), [phz(res,1);phz(res,1)+2*pi]],...
          linspace([statePlotParam{s}.hpLims,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1)',5,1));
    axis('xy');
    if s==1,                     title('Head Pitch'); end        
    if s==numel(statePlotParam), xlabel('rad');       end        
    colormap(gca,'default');
    hold('on');
    plot(statePlotParam{s}.hpLims,2*pi*ppHP.P(uind,1,s)*statePlotParam{s}.hpLims+ppHP.P(uind,2,s),'-m','LineWidth',1);
    plot(statePlotParam{s}.hpLims,2*pi*ppHP.P(uind,1,s)*statePlotParam{s}.hpLims+ppHP.P(uind,2,s)+2*pi,'-m','LineWidth',1)

    
    % head-body pitch
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(fet(res,1),[2,1]), [phz(res,1);phz(res,1)+2*pi]],...
          linspace([statePlotParam{s}.hbpLims,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1)',5,1));
    axis('xy');
    if s==1,                     title('Head-Body Pitch'); end        
    if s==numel(statePlotParam), xlabel('rad');       end        
    colormap(gca,'default');
    hold('on');
    plot(statePlotParam{s}.hbpLims,2*pi*ppHB.P(uind,1,s)*statePlotParam{s}.hbpLims+ppHB.P(uind,2,s),'-m','LineWidth',1);
    plot(statePlotParam{s}.hbpLims,2*pi*ppHB.P(uind,1,s)*statePlotParam{s}.hbpLims+ppHB.P(uind,2,s)+2*pi,'-m','LineWidth',1)
    
    % head body direction
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(circ_dist(ang(res,'head_back','head_front',1),...
                                                ang(res,'pelvis_root','spine_upper',1)),[2,1]), [phz(res,1);phz(res,1)+2*pi]],...
          linspace([statePlotParam{s}.hbdLims,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1)',5,1));
    axis('xy');
    if s==1,                     title('Head body yaw'); end        
    if s==numel(statePlotParam), xlabel('rad');       end        
    colormap(gca,'default');
    
    % drz
    subplot2(ny,nx,s,p);p=p+1;
    [out,xbins,ybins] = hist2([repmat(drz(res,unit==unitSubset),[2,1]),...
           [phz(res,1);phz(res,1)+2*pi]],...
          linspace([-1,1,50]),...
          linspace([-pi,pi*3,50]));
    imagesc(xbins,ybins,RectFilter(RectFilter(out,5,1)',5,1));
    axis('xy');    
    if s==1,                     title('DRZ');      end   
    if s==numel(statePlotParam), xlabel('AU');  end        
    colormap(gca,'default');
    hold('on');
    plot([-1,1],2*pi*ppDRZ.P(uind,1,s)*[-1,1]+ppDRZ.P(uind,2,s),'-m','LineWidth',1);
    plot([-1,1],2*pi*ppDRZ.P(uind,1,s)*[-1,1]+ppDRZ.P(uind,2,s)+2*pi,'-m','LineWidth',1)

end
drawnow();
waitforbuttonpress();
end



figure,plot([drzspk;drzspk],[phzspk;phzspk+2*pi],'.');hold('on');
plot([-1,1],2*pi*P(1,1,s)*[-1,1]+P(1,2,s),'-m','LineWidth',1);
plot([-1,1],2*pi*P(1,1,s)*[-1,1]+P(1,2,s)+2*pi,'-m','LineWidth',1)
