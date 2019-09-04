function bhv_rhm_ncp(Trial,varargin)
[mode,ncp_thresh,ncp_chan,stc_mode] = DefaultArgs(varargin,{'height',[],2,'auto_wbhr'});


disp(['bhv_rhm_ncp: ' Trial.filebase])
%% Select behavioral state collection
Trial.stc.updateMode(stc_mode);
Trial.stc.load;

%% Load Rythmic Head Motion(RHM) feature
rhm = fet_rhm(Trial);

%% Load Nasal Cavity Pressure(NCP) feature
ncp = fet_ncp(Trial,[],[],ncp_chan);

%% Whiten RHM and NCP for spectral comparison (PSD&CSD)
wang = [rhm,ncp];
wang = WhitenSignal(wang,[],1);


[ys,fs,ts] = mtcsdglong(wang,2^9,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);


figS = figure(333212);
set(figS,'Position',[268    80   685   492]);
sp(1) = subplot(311);
imagescnan({ts,fs,log10(ys(:,:,1,1))'},[-6,-3.1],0,1);axis xy,
ylabel('rhm')
sp(2) = subplot(312);
imagescnan({ts,fs,angle(ys(:,:,1,2)')},[],1,1);axis xy,
ylabel('phase diff')
sp(3) = subplot(313);
imagescnan({ts,fs,log10(ys(:,:,2,2)')},[-2,4.5],0,1);axis xy,
ylabel('ncp')
xlabel('time (s)')
linkaxes(sp,'xy');
xlim([10,40])
%Figure, 
reportfig(fullfile(Trial.path.project,'figures'),figS, ...
              'FileName',['mean_Coherence_RHM_NCP_X_' mode],...
              'Comment',[Trial.filebase ':sample:70s'])




%% Get Speed of the Body
xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(1,Trial.xyz.sampleRate));
vh = xyz.vel('spine_lower',[1,2]);


%% Construct MTADlfp to hold ys
% adjust with padding to account for spectral window
szy = size(ys);
padding = zeros([round(2^6/xyz.sampleRate/diff(ts(1:2))),szy(2:end)]);
ys = MTADlfp('data',cat(1,padding,ys,padding),'sampleRate',1/diff(ts(1:2)));


%% Resample Variable to match PSD&CSD data
vh.resample(ys);
vh.data = log10(abs(vh.data));
vh.data(~nniz(vh(:))) = nan;
xyz.resample(ys);


rhm_maxpow = MTADlfp('data',max(log10(ys(:,fs>13&fs>5,1,1)),[],2),'sampleRate',ys.sampleRate);
ncp_maxpow = MTADlfp('data',max(log10(ys(:,fs>13&fs>5,2,2)),[],2),'sampleRate',ys.sampleRate);

if isempty(ncp_thresh),
    ncp_thresh =  mean(ncp_maxpow.data);
end


    
chan_labels = {'Rhythmic Head Motion','Nasal Cavity Pressure'};

figH = figure(238482);
%set(figH,'Position',[268    80   985   692]);
set(figH,'Position',[7, 188, 1908, 1264]);
for s = 1:numel(Trial.stc.states)

    clf;
    sind = Trial.stc.states{s}.copy;
    %sind.resample(ys);
    
    switch mode
      case 'height'
        ind = nniz(xyz(sind,'head_front',3));
        ind = ncp_maxpow(sind)>ncp_thresh&ind;
        vhs = log10(xyz(sind,'head_front',3));        
        vhs = vhs(ind);
        vhlim = [1.3, 2.4];
        vh_label = 'Head Height';
        vh_units = 'mm';
      case 'speed'
        ind = nniz(vh(sind));
        ind = ncp_maxpow(sind)>ncp_thresh&ind;
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
       
    
    edges = linspace(vhlim(1),vhlim(2),15);
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


    for j = 1:size(yss,3),
        for k = 1:size(yss,3),
            subplot2(size(yss,3),size(yss,3)+1,j,k);
            if j~=k
                imagesc(1:size(edges,2),fs,vsc(:,:,j,k)'),axis xy,
                title([chan_labels{j} ' & ' chan_labels{k} ' coherence'])
            elseif j>k
                continue
            else
                imagesc(1:size(edges,2),fs,log10(vsc(:,:,j,k))'),axis xy,
                title([chan_labels{j} ' PSD Binned by ' vh_label])
            end
            set(gca,'XTickLabelMode','manual');
            set(gca,'XTickMode','manual');
            nxt = numel(get(gca,'XTickLabel'));
            nxts = cellfun(@transpose,edges_labels,'UniformOutput',false);
            xtl = strsplit(sprintf('%3.2f,',cell2mat(nxts)),',');
            xtl(end) = [];
            set(gca,'XTick',[1:2:14]);
            set(gca,'XTickLabel',xtl(1:2:end-1));
            if j~=k,caxis([.5,.9]),end
            xlabel(['Binned ' vh_label ' (' vh_units ')']);
            ylabel('Frequency Hz')
            colorbar
        end
    end


    subplot2(size(yss,3),size(yss,3)+1,2,1);
    hist2([vh(sind),log10(xyz(sind,7,3))],linspace(-1.5,1.5,30),linspace(1.3,2.3,30));
    title([Trial.filebase ' body speed Vs head height'])
    xlabel('log10 Body Speed cm/s')
    ylabel('log10 Head Height mm')
    
    nedgs = 50;
    subplot2(size(yss,3),size(yss,3)+1,1,size(yss,3)+1);
    rhm_hist_edges = prctile(rhm_maxpow(nniz(rhm_maxpow)),[5,98]);
    rhm_hist_edges = linspace(rhm_hist_edges(1),rhm_hist_edges(2),nedgs);
    hold('on')
    Na = histc(rhm_maxpow(nniz(rhm_maxpow)),rhm_hist_edges);
    nbax = bar(rhm_hist_edges,Na,'histc');axis tight,set(nbax,'FaceVertexCData',repmat([0,0,0],[nedgs,1]))
    Ns = histc(rhm_maxpow(sind),rhm_hist_edges);
    nsax = bar(rhm_hist_edges,Ns,'histc');set(nsax,'FaceVertexCData',repmat([1,.3,0],[50,1]))
    legend('all',sind.label,'Location','northeast')
    title('Max Power RHM(5-13)')
    
    subplot2(size(yss,3),size(yss,3)+1,2,size(yss,3)+1);
    ncp_hist_edges = prctile(ncp_maxpow(nniz(ncp_maxpow)),[2,98]);
    ncp_hist_edges = linspace(ncp_hist_edges(1),ncp_hist_edges(2),nedgs);
    hold('on')
    Na = histc(ncp_maxpow(nniz(ncp_maxpow)),ncp_hist_edges);
    nbax = bar(ncp_hist_edges,Na,'histc');axis tight,set(nbax,'FaceVertexCData',repmat([0,0,0],[nedgs,1]))
    Ns = histc(ncp_maxpow(sind),ncp_hist_edges);
    nsax = bar(ncp_hist_edges,Ns,'histc');set(nsax,'FaceVertexCData',repmat([1,.3,0],[50,1]))
    legend('all',sind.label,'Location','northeast')    
    Lines(ncp_thresh,[],'r');
    title('Max Power NCP(5-13)')

    reportfig(fullfile(Trial.path.project,'figures'),figH, ...
              'FileName',['mean_Coherence_RHM_NCP_X_' mode],'Comment',[Trial.filebase ':BhvState: ' ...
                        Trial.stc.states{s}.label],200  )

end

