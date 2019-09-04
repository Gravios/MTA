function [figH] = bhv_mean_coherence(Trial,varargin)
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

defargs = struct('sigRef',[],...
                 'sigCom',fet_rhm(Trial),...
                 'state','a-m',...
                 'mode',{{'height','hangle','bangle','hroll','speed','RHMpow','NCPpow'}},...
                 'ncpThresh',[],...
                 'ncpChan',2,...
                 'stcMode',Trial.stc.mode,...
                 'figHnum',238499 ...
);


[sigRef,sigCom,state,mode,ncpThresh,ncpChan,stcMode,figHnum] = DefaultArgs(varargin,defargs,'--struct');

disp(['bhv_rhm_ncp_distrb: ' Trial.filebase])

% Select behavioral state collection
stc = Trial.load('stc',stcMode);


% Load Nasal Cavity Pressure(NCP) feature if no other MTAData
% object was passed to reference signal4
if isempty(sigRef), 
    sigRef = fet_ncp(Trial,sigCom,'mta',ncpChan);
end


sig = sigRef.copy;
sig.data = [sigRef.data,...
            sigCom.data];

sparm = struct('nFFT'     ,2^9,...
               'Fs'       ,sig.sampleRate,...
               'WinLength',2^8,...
               'nOverlap' ,2^8*.875,...
               'FreqRange',[0.1,20]);

[ys,fs,ts] = fet_spec(Trial,sig,'mtcsdglong',false,[],sparm);


%% Get smoothed speed of the Body
try,
    xyz = Trial.load('xyz','seh');
catch
    xyz = Trial.load('xyz');
end

xyz.filter('ButFilter',3,2.4,'low');

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
ncp_maxpow = MTADlfp('data',max(nys(:,fs>13&fs>5,2,2),[],2),'sampleRate',ys.sampleRate);

    
chan_labels = {sigRef.label,sigCom.label};

set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[0,0,(4+2)*numel(mode)+2,4*(3+2)+2])
set(hfig,'PaperPositionMode','auto');
for s = 1;%:numel(stc.states)

    clf;
    sind = stc{state,xyz.sampleRate};
    sind.data(sind.data(:,1)>xyz.size(1),:)=[];
    sind.data(sind.data(:,2)>xyz.size(1),:)=[];

    for m = 1:numel(mode),
        switch mode{m}
          case 'height'
            ind = nniz(xyz(sind,'head_front',3));
            %ind = ncp_maxpow(sind)>ncpThresh&ind;
            vhs = log10(xyz(sind,'head_front',3));        
            vhs = vhs(ind);
            vhlim = [0, 2.5];
            vh_label = 'Height(H)';
            vh_units = 'log10(mm)';
          case 'hangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'head_back','head_front',2));
            %ind = ncp_maxpow(sind)>ncpThresh&ind;
            vhs = ang(sind,'head_back','head_front',2);
            vhs = vhs(ind);
            vhlim = [-1.5, 1.5];
            vh_label = 'Pitch(H)';
            vh_units = 'rad';
          case 'bangle'
            ang = Trial.ang.copy;
            ang = ang.create(Trial,xyz);
            ind = nniz(ang(sind,'spine_lower','spine_upper',2));
            %ind = ncp_maxpow(sind)>ncpThresh&ind;
            vhs = ang(sind,'spine_lower','spine_upper',2);
            vhs = vhs(ind);
            vhlim = [0, .6];
            vh_label = 'Pitch(B)';
            vh_units = 'rad';
          case 'hroll'
            hbflr = transform_origin(Trial,xyz,'head_back','head_front',{'head_left','head_right'});            
            hroll = MTADxyz('data',hbflr.roll,'sampleRate',xyz.sampleRate);
            ind = nniz(ang(sind,'spine_lower','spine_upper',2));
            %ind = ncp_maxpow(sind)>ncpThresh&ind;
            vhs = hroll(sind);
            vhs = vhs(ind);
            vhlim = [-pi/2,pi/2];
            vh_label = 'roll(H)';
            vh_units = 'rad';
          case 'speed'
            ind = nniz(vh(sind));
            %ind = ncp_maxpow(sind)>ncpThresh&ind;
            vhs =vh(sind);
            vhs = vhs(ind);
            %vhlim =prctile(vhs,[5,98]);
            vhlim = [-3,2.2];
            vh_label = 'Speed(B)';
            vh_units = 'log10(cm/s)';
          case 'NCPpow'
            vhncp = MTADlfp('data',ncp_maxpow(sind),'sampleRate',ncp_maxpow.sampleRate);
            ind = nniz(vh(sind));
            vhs =vhncp;
            vhs = vhs(ind);
            vhlim =prctile(vhs,[2,98]);
            vh_label = 'NCP(5-13)';      
            vh_units = 'log10(mV^2/s)';
          case 'RHMpow'
            vhrhm = rhm_maxpow(sind);
            ind = nniz(vh(sind));
            vhs = clip(vhrhm,-9,5);
            vhs = vhs(ind);
            vhlim =prctile(vhs,[2,98]);
            vh_label = 'RHM(5-13)';
            vh_units = 'log10(cm^2/s)';
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

        vsc = nan([size(edges,2),numel(fs),size(yss,3),size(yss,3)]);
        for i = edges,
            ind = i(1)<=vhs&vhs<i(2);
            for j = 1:size(yss,3),
                for k = 1:size(yss,3),
                    if sum(ind)~=0,
                        if k~=j,
                            vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))...
                                                                ./nanmean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
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
        axes('Units','centimeters',...
             'Position',[(4+2)*(m-1)+2,(3+2)*4-(4)+1,4,3])
        
        %imagesc(vedgs,fs,mrv(:,:,1)',[.1,max(prctile(mrv(nniz(mrv),:,1),[95]),[],2)]),axis xy      
        
        imagesc(vedgs,fs,bsxfun(@plus,mrv(:,:,1)',abs(min(mrv(:,:,1)'))),...
                [.1,max(prctile(mrv(nniz(mrv),:,1),[95]),[],2)]),axis xy
        title([chan_labels{1} ' mean PSD given ' vh_label])
        xlabel([ vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz')
        colormap jet
        caxis(prctile(reshape(get(get(gca,'Children'),'CData'),[],1),[5,95]))
        h=colorbar('EastOutside');
        apos = get(gca,'Position');
        set(h,'Units','centimeters');
        set(h,'Position',[apos(1)+4.1,apos(2),0.3,3])

        %% NCP psd
        axes('Units','centimeters',...
             'Position',[(4+2)*(m-1)+2,(3+2)*3-(4)+1,4,3])
        imagesc(vedgs,fs,bsxfun(@plus,mrv(:,:,2)',abs(min(mrv(:,:,2)'))),...
                [.1,max(prctile(mrv(nniz(mrv),:,2),[95]),[],2)]),axis xy
        %imagesc(vedgs,fs,mrv(:,:,2)',[.1,max(prctile(mrv(nniz(mrv),:,2),[95]),[],2)]);
        axis('xy');
        title([chan_labels{2} ' mean PSD given ' vh_label]);
        xlabel([ vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz');
        colormap('jet');
        caxis(prctile(reshape(get(get(gca,'Children'),'CData'),[],1),[5,95]));
        h=colorbar('EastOutside');
        apos = get(gca,'Position');
        set(h,'Units','centimeters');
        set(h,'Position',[apos(1)+4.1,apos(2),0.3,3]);


        %% RHM NCP coherence
        axes('Units','centimeters',...
             'Position',[(4+2)*(m-1)+2,(3+2)*2-(4)+1,4,3])
        imagesc(edges(1,:),fs,vsc(:,:,1,2)');        
        axis('xy');
        title([chan_labels{1} ' & ' chan_labels{2} ' coherence']);
        xlabel([ vh_label ' (' vh_units ')']);
        ylabel('Frequency Hz');
        caxis([0.4,1]);
        h=colorbar('EastOutside');
        apos = get(gca,'Position');
        set(h,'Units','centimeters');
        set(h,'Position',[apos(1)+4.1,apos(2),0.3,3])

        
        %subplot2(4,numel(mode),4,m);
        axes('Units','centimeters',...
             'Position',[(4+2)*(m-1)+2,(3+2)*1-(4)+1,4,3])
        bar(vedgs,N,'histc')
        ylabel('Count')
        xlabel(['Binned ' vh_label ' (' vh_units ')']);
        xlim(edges([1,end]))

    end

    figPath = fullfile(Trial.spath,'figures');
    if ~exist(figPath,'dir'),
        mkdir(figPath);
    end
    
    saveas(hfig, fullfile(figPath,[mfilename '.fig']),'fig');


% $$$     reportfig(Trial.spath,    ... Path where figures are stored
% $$$               hfig,           ... Figure handle
% $$$               [mfilename],    ... Figure Set Name
% $$$               'figures',      ... Directory where figures reside
% $$$               false,          ... Do Not Preview
% $$$               '',             ... thmb_cap  
% $$$               '',             ... exp_cap
% $$$               [],             ... Resolution
% $$$               true,           ... Save FIG
% $$$               'png',(4+2)*numel(mode)+2,4*(3+2)+2);       % Output Format
    

% $$$     reportfig(fullfile(Trial.path.project,'figures'),hfig, ...
% $$$               'FileName',['mean_Coherence_RHM_NCP_X_' strjoin(mode,'_')],...
% $$$               'Comment',[Trial.filebase ':BhvState: ' stc.states{s}.label],...
% $$$               'Resolution',200,...
% $$$               'SaveFig',true)
end