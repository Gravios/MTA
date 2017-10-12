function [figH] = bhv_rhm_ncp_distrb(Trial,varargin)
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

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('mode',           {{,'hangle','bhangle','bspeed','hspeed','NCPpow','RHMpow'}},  ...
                 'ncpThreshold',   [],                                                           ...
                 'ncpChannel',     2,                                                            ...
                 'stcMode',        'msnn_ppsvd',                                                 ...
                 'p',              8,                                                            ...
                 'xyzProcOpts',    {'SPLINE_SPINE_HEAD_EQD'}                                     ...
);

[mode,ncpThreshold,ncpChannel,stcMode,p,xyzProcOpts] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------




xyz = preproc_xyz(Trial,xyzProcOpts);
xyz.filter('ButFilter',3,2.4,'low');

disp(['bhv_rhm_ncp_distrb: ' Trial.filebase])

% SELECT behavioral state collection
stc = Trial.load('stc',stcMode);


% LOAD Rythmic Head Motion(RHM) feature
%rhm = fet_rhm_exp(Trial);
%rhm = fet_rbm(Trial);
%rhm = fet_rhm(Trial);
%rhm = fet_rhmPCA(Trial);
rhm = fet_rhm(Trial);
%rhm = fet_rhp(Trial);
%rhm = fet_rhm_rf(Trial);
rhm.resample(xyz);
rhm.data(~nniz(rhm.data(:))) = 1;

% LOAD Nasal Cavity Pressure(NCP) feature
ncp = fet_ncp(Trial,rhm,'mta',ncpChannel);


% WHITEN RHM and NCP for spectral comparison (PSD&CSD)
% $$$ wang = [rhm.data,ncp.data];
% $$$ wang = WhitenSignal(wang,[],1);

rhm.data = [rhm.data,ncp.data];


sparm = struct('nFFT'     ,2^(p),...
               'Fs'       ,rhm.sampleRate,...
               'WinLength',2^(p-1),...
               'nOverlap' ,2^(p-1)*.875,...
               'FreqRange',[.5,20]);

[ys,fs,ts] = fet_spec(Trial,rhm,'mtcsdglong',false,[],sparm);


% GET smoothed speed of the Body
vh = xyz.vel({'spine_lower','head_front'},[1,2]);
vh.resample(ys);
vh.data = log10(abs(vh.data));
vh.data(~nniz(vh(:))) = nan;

xyz.resample(ys);

nys = ys.copy;
nys.data = log10(nys.data);
nys.data(nys<-9) = nan;
%nys.data(nys>5) = nan;
nys.data(~nniz(nys.data))=nan;
nys.data = (nys.data-repmat(nanmedian(nys(nniz(nys),:,:,:)),[nys.size(1),1,1]))./repmat(nanstd(nys(nniz(nys),:,:,:)),[nys.size(1),1,1]);

rhm_maxpow = MTADlfp('data',max(nys(:,fs>13&fs>5,1,1),[],2),'sampleRate',ys.sampleRate);
ncp_maxpow = MTADlfp('data',max(nys(:,fs>13&fs>5,2,2),[],2),'sampleRate',ys.sampleRate);

    
chan_labels = {'RHM','NCP'};

figH = figure(gen_figure_id);
figH.Units = 'centimeters';
figH.Position = [1,1,30,20];
for s = 1;

    clf;
    %sind = stc.states{s}.copy;
    sind = stc{'a-m-s'};%+[1,-1];
    %sind = stc{'w+n+p'};%+[1,-1];
    %sind = stc{'a&w'};%+[1,-1];
    %sind = stc{'w+q'}+[1,-1];
    %sind.resample(ys);

    for m = 1:numel(mode),
        switch mode{m}
          case 'height'
            ind = nniz(xyz(sind,'head_front',3));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs = log10(xyz(sind,'head_front',3));
            vhs = vhs(ind);
            vhlim = [1, 2.2];
            vh_label = 'Head Height';
            vh_units = 'mm';
          case 'hangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'head_back','head_front',2));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs = ang(sind,'head_back','head_front',2);
            vhs = vhs(ind);
            vhlim = [-1.5, 1.5];
            vh_label = 'Head Pitch';
            vh_units = 'rad';
          case 'bangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'spine_lower','spine_upper',2));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs = ang(sind,'spine_lower','spine_upper',2);
            vhs = vhs(ind);
            vhlim = [-1.5, 1.5];
            vh_label = 'Body Pitch';
            vh_units = 'rad';
          case 'bhangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'spine_upper','head_front',2));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs = ang(sind,'spine_lower','spine_upper',2);
            vhs = vhs(ind);
            vhlim = [-1.5, 1.5];
            vh_label = 'Body Pitch';
            vh_units = 'rad';
          case 'bspeed'
            ind = nniz(vh(sind));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs =vh(sind,1);
            vhs = vhs(ind);
            %vhlim =prctile(vhs,[5,98]);
            vhlim = [-3,2];
            vh_label = 'Body Speed';
            vh_units = 'cm/s';
          case 'hspeed'
            ind = nniz(vh(sind));
            %ind = ncp_maxpow(sind)>ncpThreshold&ind;
            vhs =vh(sind,2);
            vhs = vhs(ind);
            %vhlim =prctile(vhs,[5,98]);
            vhlim = [-3,2];
            vh_label = 'Head Speed';
            vh_units = 'cm/s';
          case 'NCPpow'
            vhncp = MTADlfp('data',ncp_maxpow(sind),'sampleRate',ncp_maxpow.sampleRate);
            ind = nniz(vh(sind));
            vhs =vhncp;
            vhs = vhs(ind);
            vhlim =prctile(vhs,[2,98]);
            vh_label = 'NCP_pow(5-13)';      
            vh_units = 'mV^2/s';
          case 'RHMpow'
            vhrhm = rhm_maxpow(sind);
            ind = nniz(vh(sind));
            vhs = clip(vhrhm,-9,5);
            vhs = vhs(ind);
            vhlim =prctile(vhs,[2,98]);
            vh_label = 'RHM_pow(5-13)';
            vh_units = 'cm^2/s';
        end 
        vys = nys(sind,:,:,:);
        vys = vys(ind,:,:,:);


        vbins = 30;
        %vedgs = linspace(vhlim(1),vhlim(2),vbins);
        vedgs = prctile(vhs,linspace(1,99,vbins));
        [N,vbs] = histc(vhs,vedgs);
        mrv = nan(numel(N),nys.size(2),nys.size(3));
        for c = 1:nys.size(3),
            for f =1:nys.size(2),
                mrv(:,f,c) = accumarray(vbs(nniz(vbs)&nniz(vys)),...
                                        vys(nniz(vbs)&nniz(vys),f,c,c),...
                                        [vbins,1],...
                                        @nanmean);
            end
        end

        
        edges = prctile(vhs,linspace(1,99,vbins));        
        %edges = linspace(vhlim(1),vhlim(2),30);
        edges = [edges(1:end-1);edges(2:end)];
        edges_labels = mat2cell(10.^mean(edges),1,ones(1,size(edges,2)))';
        yss = ys(sind,:,:,:);
        yss = yss(ind,:,:,:);
        clear ind;

        vsc = [];
        for i = edges,
            ind = i(1)>=vhs&vhs<i(2);
            if sum(ind)<10, continue,end
            for j = 1:size(yss,3),
                for k = 1:size(yss,3),
                    if sum(ind)~=0,
                        if k~=j,
                            vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))./nanmean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
                        else
                            vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(yss(ind,:,k,j));
                        end

                    else
                        vsc(find(i(1)==edges(1,:)),:,k,j) = zeros([1,size(yss,2)]);
                    end

                end
            end
        end

        
        %% RHM psd
        subplot2(3,numel(mode),1,m);
        imagesc(vedgs,fs,mrv(:,:,1)',[0,max(prctile(mrv(nniz(mrv),:,1),[95]),[],2)]);,axis xy
        title({[chan_labels{1} ' mean PSD Binned by '],[ vh_label]})
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        h = colorbar;
        h.Position(1) = h.Position(1) +.05;
        caxis([-1.5,1.5])

        %% NCP psd
        subplot2(3,numel(mode),2,m);        
        imagesc(vedgs,fs,mrv(:,:,2)',[0,max(prctile(mrv(nniz(mrv),:,2),[95]),[],2)]);,axis xy
        title({[chan_labels{2} ' mean PSD Binned by '],[ vh_label]})
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        h = colorbar;
        h.Position(1) = h.Position(1) +.05;
        caxis([-1.5,1.5])

        %% RHM NCP coherence
        subplot2(3,numel(mode),3,m);
        imagesc(edges(1,:),fs,vsc(:,:,1,2)');,axis xy,
        title({[chan_labels{1} ' & '],[ chan_labels{2} ' coherence']})
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        h = colorbar;
        h.Position(1) = h.Position(1) +.05;
        caxis([0.5,1]);

% $$$         subplot2(numel(mode),4,m,4);
% $$$         bar(vedgs,N,'histc')
% $$$         ylabel('Count')
% $$$         xlabel(['Binned ' vh_label ' (' vh_units ')']);

    end
    colormap jet;

% $$$     reportfig(fullfile(Trial.path.data,'figures'),figH, ...
% $$$               'FileName',['mean_Coherence_RHM_NCP_X_' strjoin(mode,'_')],...
% $$$               'Comment',[Trial.filebase ':BhvState: ' stc.states{s}.label],...
% $$$               'Resolution',200,...
% $$$               'SaveFig',true)
end