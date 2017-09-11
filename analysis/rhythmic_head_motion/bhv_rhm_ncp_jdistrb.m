function [figH] = bhv_rhm_ncp_jdistrb(varargin)
%function bhv_rhm_ncp_distrb(Trial,varargin)
%
%
%  varargin:
%
%    mode:    string/cellArray - predetermined 
%             Def({'height','hangle'})
%
%
%    ncp_thresh: duno    - Def([]), duno what it's for
%
%    ncp_chan:   numeric - Def(2), 
%
%    stcMode:   string  - Def('auto_wbhr'), 



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',    'ncp',                                                        ...
                 'ncpChannel',     2,                                                            ...
                 'stcMode',        'msnn_ppsvd',                                                 ...
                 'p',              8,                                                            ...
                 'xyzProcOpts',    {{'SPLINE_SPINE_HEAD_EQI'}}                                   ...
);

[sessionList,ncpChannel,stcMode,p,xyzProcOpts] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% GET and CONVERT structArray list to cellArray of structs
slist = subsref(cf(@(x)  mat2cell(x,1,ones(size(x))),  {get_session_list(sessionList)}),substruct('{}',{1}));


if iscell(sessionList) & isa(sessionList{1},'MTASession'),
% ASSIGN cell array of MTATrials for behavioral labeling    
    Trials = sessionList;
elseif isa(sessionList,'MTATrial'),
% ENCAPSULATE trial in cell
    Trials = {sessionList};
else
% LOAD Trials for behavioral labeling
    Trials = cf(@(t)  MTATrial.validate(t), slist);
end
numTrials = numel(slist);



% LOAD xyz
xyz = cf(@(t,o)  preproc_xyz(t,o)  ,Trials,repmat({xyzProcOpts},[1,numTrials]));
      cf(@(x)    x.filter('ButFilter',3,2.4,'low'),  xyz);      

% LOAD the state collections
if ~isempty(stcMode) && ischar(stcMode),
    stc = cf(@(t,m)  t.load('stc',m),  Trials,repmat({stcMode},[1,numTrials]));
else
    stc = cf(@(t)    t.load('stc'),    Trials);
end



% LOAD Rythmic Head Motion(RHM) feature
rhm = cf(@(t)    fet_rhm(t),     Trials);
      cf(@(r,x)  r.resample(x),  rhm, xyz);
for t = 1:numTrials, rhm{t}.data(~nniz(rhm{t}.data(:))) = 1; end

% LOAD Nasal Cavity Pressure(NCP) feature
ncp = cf(@(t,r,s)  fet_ncp(t,r,'mta',s.ncpChannel),  Trials, rhm, slist);
ncp9hz = cf(@(n) n.copy ncp

% CREATE fe
fet =  cf(@(r)      r.copy(),                       rhm          );
       cf(@(f,r,n)  set(f,'data',[r.data,n.data]),  fet, rhm, ncp);

% SET spectral estimation parameters
sparm = struct('nFFT'     ,2^(p),...
               'Fs'       ,fet{1}.sampleRate,...
               'WinLength',2^(p-1),...
               'nOverlap' ,2^(p-1)*.875,...
               'FreqRange',[.5,20]);

% COMPUTE auto and cross spectrum
[ys,fs,ts] = cf(@(t,f,p) fet_spec(t,f,'mtcsdglong',false,[],p),  Trials,fet,repmat({sparm},[1,numTrials]));


% GET smoothed speed of the Body
fxyz = cf(@(x,y)  resample(copy(x),y),                  xyz, ys);
vh   = cf(@(x)    x.vel({'spine_lower','hcom'},[1,2]),  fxyz   );
       cf(@(v)    set(v,'data',log10(abs(v.data))),     vh     );
for t = 1:numTrials, vh{t}.data(~nniz(vh{t}.data(:))) = nan; end

    
% LOAD head body pitch features
features = cf(@(t)    fet_HB_pitch(t),  Trials      );
           cf(@(f,y)  f.resample(y),    features, ys);

% LOAD periods of state subset
sper = cf(@(s)   [s{'w+n+p+r&a'}],      stc          );
       cf(@(p)   p.cast('TimeSeries'),  sper         );
       cf(@(p)   set(p,'data',logical(p.data)), sper );
       cf(@(p,y) p.resample(y),         sper, ys     );


for s = 1:numTrials,       

    afet{s} = features{s}.data;
    aper{s} = sper{s}.data;
    nbins = 20;


    fetInds{s} = [];
    eds{s} = [];
    for f = 1:size(features{1},2),
        fetSub{s} = afet{s}(aper{s},f);    
        %eds(:,f) = prctile(fetSub(randi([1,sum(sper.data)],[1e6,1]),1),linspace(0.1,99.9,nbins));
        eds{s}(:,f) = linspace([prctile(fetSub{s}(randi([1,sum(aper{s})],[1e6,1]),1),[0.1,99.9]),nbins]);
        %eds(:,f) = prctile(fetSub(randi([1,sum(sper.data)],[1e6,1]),1),linspace(0.1,99.9,nbins));
        fetInds{s}(:,f) = discretize(fetSub{s},eds{s}(:,f));
    end

end


% spectral analysis over sessions and frequency band
for s = 1:numTrials,       

    ays{s}   = ys{s}.data;
    ysegs{s} = mean(ays{s}(aper{s},fs{1}<12&fs{1}>6,:,:),2);

    fetInds{s}(~nniz(fetInds{s}),:) = [];
    ysegs{s}(~nniz(fetInds{s}),:,:,:) = [];


    rout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
    pout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
    mout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
    sout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
    mouts{s} = nan([nbins-1,nbins-1,numel(fs{1})]);
    souts{s} = nan([nbins-1,nbins-1,numel(fs{1})]);
    out{s}   = nan([nbins-1,nbins-1,numel(fs{1})]);
    freqInd = 1;
    
    for f = 1:nbins-1,
        for g = 1:nbins-1,
            ind = fetInds{s}(:,1)==f & fetInds{s}(:,2)==g;
            if sum(ind)>5,
                out{s}(f,g,freqInd) = nanmean(abs(ysegs{s}(ind,freqInd,1,2)))./...
                    nanmean(sqrt(ysegs{s}(ind,freqInd,1,1).*ysegs{s}(ind,freqInd,2,2)));
                pout{s}(f,g,freqInd) = circ_mean(angle(ysegs{s}(ind,freqInd,1,2)));
                rout{s}(f,g,freqInd) = abs(mean(ysegs{s}(ind,freqInd,1,2)));
                mout{s}(f,g,freqInd) = mean(abs(ysegs{s}(ind,freqInd,1,1)));
                sout{s}(f,g,freqInd) = mean(abs(ysegs{s}(ind,freqInd,2,2)));            
                mouts{s}(f,g,freqInd) = std(abs(ysegs{s}(ind,freqInd,1,1)));
                souts{s}(f,g,freqInd) = std(abs(ysegs{s}(ind,freqInd,2,2)));            
                %pout(f,g,freqInd) = angle(mean(imag(ysegs(ind,freqInd,1,2))));
            end    
        end
    end

end



s = 4;
ieds = diff(eds{s})+eds{s}(1:end-1,:);
figure();
for s = 1:numTrials,    
    subplot(3,4,s);
    imagescnan({ieds(:,1),ieds(:,2),out{s}(:,:,1)'},[0,1],'linear');
    axis('xy');
    title(['coherence: ',num2str([6,12])]);    
end

figure();
for s = 1:numTrials,    
    subplot(3,4,s);
    imagescnan({ieds(:,1),ieds(:,2),pout{s}(:,:,1)'},[-pi,pi],'circular',false,[],1,1,@hsv);
    axis('xy');
    title(['phase diff: ',num2str([6,12])]);
end





% spectral analysis over sessions and frequency bins
for s = 1:numTrials,       


    ays{s}   = ys{s}.data;
    ysegs{s} = ays{s}(aper{s},:,:,:);

    fetInds{s}(~nniz(fetInds{s}),:) = [];
    ysegs{s}(~nniz(fetInds{s}),:,:,:) = [];

    [freq,freqInd] = NearestNeighbour(fs{1},7);

rout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
pout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
mout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
sout{s}  = nan([nbins-1,nbins-1,numel(fs{1})]);
mouts{s} = nan([nbins-1,nbins-1,numel(fs{1})]);
souts{s} = nan([nbins-1,nbins-1,numel(fs{1})]);
out{s}   = nan([nbins-1,nbins-1,numel(fs{1})]);
for freqInd = 1:numel(fs{1}),
    for f = 1:nbins-1,
        for g = 1:nbins-1,
            ind = fetInds{s}(:,1)==f & fetInds{s}(:,2)==g;
            if sum(ind)>5,
                out{s}(f,g,freqInd) = nanmean(abs(ysegs{s}(ind,freqInd,1,2)))./...
                    nanmean(sqrt(ysegs{s}(ind,freqInd,1,1).*ysegs{s}(ind,freqInd,2,2)));
                pout{s}(f,g,freqInd) = circ_mean(angle(ysegs{s}(ind,freqInd,1,2)));
                rout{s}(f,g,freqInd) = abs(mean(ysegs{s}(ind,freqInd,1,2)));
                mout{s}(f,g,freqInd) = mean(abs(ysegs{s}(ind,freqInd,1,1)));
                sout{s}(f,g,freqInd) = mean(abs(ysegs{s}(ind,freqInd,2,2)));            
                mouts{s}(f,g,freqInd) = std(abs(ysegs{s}(ind,freqInd,1,1)));
                souts{s}(f,g,freqInd) = std(abs(ysegs{s}(ind,freqInd,2,2)));            
                %pout(f,g,freqInd) = angle(mean(imag(ysegs(ind,freqInd,1,2))));
            end    
        end
    end
end

end

s = 4;
ieds = diff(eds{s})+eds{s}(1:end-1,:);
figure();
indexF = 6:3:24,

for s = 1:numTrials,    
for f = indexF,    
    subplot2(numTrials,numel(indexF),s,find(f==indexF));
    imagescnan({ieds(:,1),ieds(:,2),out{s}(:,:,f)'},[0,1],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end
end

figure();
for s = 1:numTrials,
for f= indexF,
    subplot2(numTrials,numel(indexF),s,find(f==indexF));    
    imagescnan({ieds(:,1),ieds(:,2),pout{s}(:,:,f)'},[-pi,pi],'circular',false,[],1,1,@hsv);
    axis('xy');
    title(num2str(fs{1}(f)));
end
end

% $$$ ays = cf(@(y)  y.data,  ys);
% $$$ ays = cat(1, ays{:});
% $$$ 
% $$$ ysegs = ays(aper,:,:,:);
% $$$ 
% $$$ fetInds(~nniz(fetInds),:) = [];
% $$$ ysegs(~nniz(fetInds),:,:,:) = [];
% $$$ 
% $$$ [freq,freqInd] = NearestNeighbour(fs{1},7);
% $$$ 
% $$$ rout  = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ pout  = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ mout  = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ sout  = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ mouts = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ souts = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ out   = nan([nbins-1,nbins-1,numel(fs{1})]);
% $$$ for freqInd = 1:numel(fs{1}),
% $$$     for f = 1:nbins-1,
% $$$         for g = 1:nbins-1,
% $$$             ind = fetInds(:,1)==f & fetInds(:,2)==g;
% $$$             if sum(ind)>5,
% $$$                 out(f,g,freqInd) = nanmean(abs(ysegs(ind,freqInd,1,2)))./...
% $$$                     nanmean(sqrt(ysegs(ind,freqInd,1,1).*ysegs(ind,freqInd,2,2)));
% $$$                 pout(f,g,freqInd) = circ_mean(angle(ysegs(ind,freqInd,1,2)));
% $$$                 rout(f,g,freqInd) = abs(mean(ysegs(ind,freqInd,1,2)));
% $$$                 mout(f,g,freqInd) = mean(abs(ysegs(ind,freqInd,1,1)));
% $$$                 sout(f,g,freqInd) = mean(abs(ysegs(ind,freqInd,2,2)));            
% $$$                 mouts(f,g,freqInd) = std(abs(ysegs(ind,freqInd,1,1)));
% $$$                 souts(f,g,freqInd) = std(abs(ysegs(ind,freqInd,2,2)));            
% $$$                 %pout(f,g,freqInd) = angle(mean(imag(ysegs(ind,freqInd,1,2))));
% $$$             end    
% $$$         end
% $$$     end
% $$$ end




% $$$ figure,
% $$$ imagesc(eds(:,1),eds(:,2),out)
% $$$ axis xy

ieds = diff(eds)+eds(1:end-1,:);

figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),out(:,:,f)'},[0,1],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end

figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),pout(:,:,f)'},[-pi,pi],'circular',false,[],1,1,@hsv);
    axis('xy');
    title(num2str(fs{1}(f)));
end


figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),log10(rout(:,:,f))'},[-2,2],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end

figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),log10(mout(:,:,f))'},[-4,-2],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end


figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),mouts(:,:,f)'},[],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end


figure();
for f= 1:30
    subplot(5,6,f);
    imagescnan({ieds(:,1),ieds(:,2),log10(sout(:,:,f))'},[3,7],'linear');
    axis('xy');
    title(num2str(fs{1}(f)));
end





    

    
% $$$ 
% $$$     figH = figure(48849);
% $$$     plot(anq.(nqFields{1})(anq.eDist>eDistThreshold),anq.(nqFields{2})(anq.eDist>eDistThreshold),'.')
% $$$     title('Fit for hyperplane perpendicular to selected dimensions')
% $$$     xlabel(nqFields{1});
% $$$     ylabel(nqFields{2});
% $$$     pram = draw_lines(figH,'line_fit');
% $$$     %saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.png');
% $$$     %saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.fig');
% $$$ 
% $$$     usp.pram   = pram;
% $$$     usp.fields = nqFields;
% $$$ 
% $$$     save(fullfile(Trial.path.cfg,'unit_selection_criteria.mat'),'usp','-v7.3');
% $$$     units = [];
