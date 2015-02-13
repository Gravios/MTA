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
%    stc_mode:   string  - Def('auto_wbhr'), 

[mode,ncp_thresh,ncp_chan,stc_mode] = DefaultArgs(varargin,{{'height','hangle'},[],2,'auto_wbhr'});


disp(['bhv_rhm_ncp_distrb: ' Trial.filebase])
%% Select behavioral state collection
Trial.stc.updateMode(stc_mode);
Trial.stc.load;

%% Load Rythmic Head Motion(RHM) feature
rhm = fet_rhm(Trial,[],'default');

%% Load Nasal Cavity Pressure(NCP) feature
ncp = fet_ncp(Trial,rhm,'default',ncp_chan);

%% Whiten RHM and NCP for spectral comparison (PSD&CSD)
% $$$ wang = [rhm.data,ncp.data];
% $$$ wang = WhitenSignal(wang,[],1);

rhm.data = [rhm.data,ncp.data];

sparm = struct('nFFT'     ,2^9,...
               'Fs'       ,rhm.sampleRate,...
               'WinLength',2^7,...
               'nOverlap' ,2^7*.875,...
               'FreqRange',[1,20]);

[ys,fs,ts] = fet_spec(Trial,rhm,'mtcsdglong',false);


%% Get smoothed speed of the Body
xyz = Trial.load('xyz').filter(gtwin(1,Trial.xyz.sampleRate));

vh = xyz.vel('spine_lower',[1,2]);
vh.resample(ys);
vh.data = log10(abs(vh.data));
vh.data(~nniz(vh(:))) = nan;

xyz.resample(ys);

nys = ys.copy;
nys.data = log10(nys.data);
nys.data(nys<-9) = nan;
nys.data(~nniz(nys.data))=nan;
nys.data = (nys.data-repmat(nanmedian(nys(nniz(nys),:,:,:)),[nys.size(1),1,1]))./repmat(nanstd(nys(nniz(nys),:,:,:)),[nys.size(1),1,1]);

rhm_maxpow = MTADlfp('data',max(nys(:,fs>13&fs>5,1,1),[],2),'sampleRate',ys.sampleRate);
ncp_maxpow = MTADlfp('data',max(nys(:,fs>12&fs>6,2,2),[],2),'sampleRate',ys.sampleRate);

    
chan_labels = {'Rhythmic Head Motion','Nasal Cavity Pressure'};

figH = figure(238482);
%set(figH,'Position',[268    80   985   692]);
pos =[67,441,1607,482];
%pos =[20, 452, 1642, 471];
set(figH,'Position',pos);
for s = 1;%:numel(Trial.stc.states)

    clf;
    sind = Trial.stc.states{s}.copy;
    %sind.resample(ys);

    for m = 1:numel(mode),
        switch mode{m}
          case 'height'
            ind = nniz(xyz(sind,'head_front',3));
            %ind = ncp_maxpow(sind)>ncp_thresh&ind;
            vhs = log10(xyz(sind,'head_front',3));        
            vhs = vhs(ind);
            vhlim = [1, 2.2];
            vh_label = 'Head Height';
            vh_units = 'mm';
          case 'hangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'head_back','head_front',2));
            %ind = ncp_maxpow(sind)>ncp_thresh&ind;
            vhs = ang(sind,'head_back','head_front',2);
            vhs = vhs(ind);
            vhlim = [-1.5, 1.5];
            vh_label = 'Head Angle';
            vh_units = 'rad';
          case 'speed'
            ind = nniz(vh(sind));
            %ind = ncp_maxpow(sind)>ncp_thresh&ind;
            vhs =vh(sind);
            vhs = vhs(ind);
            %vhlim =prctile(vhs,[5,98]);
            vhlim = [-1.5,1.5];
            vh_label = 'Body Speed';
            vh_units = 'cm/s';
          case 'NCPpow'
            vhncp = MTADlfp('data',ncp_maxpow(sind),'sampleRate',ncp_maxpow.sampleRate);
            ind = nniz(vh(sind));
            vhs =vhncp;
            vhs = vhs(ind);
            vhlim =prctile(vhs,[5,95]);
            vh_label = 'NCP_pow(5-13)';      
            vh_units = 'mV^2/s';
          case 'RHMpow'
            vhrhm = rhm_maxpow(sind);
            ind = nniz(vh(sind));
            vhs = clip(vhrhm,-9,5);
            vhs = vhs(ind);
            vhlim =prctile(vhs,[5,95]);
            vh_label = 'RHM_pow(5-13)'
            vh_units = 'cm^2/s'
        end 
        vys = nys(sind,:,:,:);
        vys = vys(ind,:,:,:);


        vbins = 30;
        vedgs = linspace(vhlim(1),vhlim(2),vbins);
        [N,vbs] = histc(vhs,vedgs);
        mrv = nan(numel(N),nys.size(2),nys.size(3));
        for c = 1:nys.size(3),
        for f =1:nys.size(2),
            mrv(:,f,c) = accumarray(vbs(nniz(vbs)&nniz(vys)),vys(nniz(vbs)&nniz(vys),f,c,c),[vbins,1],@nanmean);
        end
        end

        
        edges = linspace(vhlim(1),vhlim(2),30);
        edges = [edges(1:end-1);edges(2:end)];
        edges_labels = mat2cell(10.^mean(edges),1,ones(1,size(edges,2)))';
        yss = ys(sind,:,:,:);
        yss = yss(ind,:,:,:);
        clear ind;

        vsc = [];
        for i = edges,
            ind = i(1)>=vhs&vhs<i(2);
            for j = 1:size(yss,3),
                for k = 1:size(yss,3),
                    if sum(ind)~=0,
                        if k~=j,
                            vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))./mean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
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
        subplot2(numel(mode),4,m,1);
        imagesc(vedgs,fs,mrv(:,:,1)',[.1,max(prctile(mrv(nniz(mrv),:,1),[95]),[],2)]),axis xy
        title([chan_labels{1} ' mean PSD Binned by ' vh_label])
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        colorbar

        %% NCP psd
        subplot2(numel(mode),4,m,2);
        imagesc(vedgs,fs,mrv(:,:,2)',[.1,max(prctile(mrv(nniz(mrv),:,2),[95]),[],2)]),axis xy
        title([chan_labels{2} ' mean PSD Binned by ' vh_label])
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        colorbar

        %% RHM NCP coherence
        subplot2(numel(mode),4,m,3);
        imagesc(edges(1,:),fs,vsc(:,:,1,2)'),axis xy,
        title([chan_labels{1} ' & ' chan_labels{2} ' coherence'])
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        colorbar

        
        subplot2(numel(mode),4,m,4);
        bar(vedgs,N,'histc')
        ylabel('Count')
        xlabel(['Binned ' vh_label ' (' vh_units ')']);

    end

    reportfig(fullfile(Trial.path.data,'figures'),figH, ...
              'FileName',['mean_Coherence_RHM_NCP_X_' strjoin(mode,'_')],...
              'Comment',[Trial.filebase ':BhvState: ' Trial.stc.states{s}.label],...
              'Resolution',200,...
              'SaveFig',true)
end