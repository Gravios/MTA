function nq = NeuronQuality(Session, varargin)
%function NeuronQuality(Session, Electrodes,Display,Batch,overwrite)
%
%  Compute statistics on the qualities of neuronal unit clusters for a session
%  
%

Par = LoadXml(fullfile(Session.spath,[Session.name '.xml']));

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Electrodes',                      [1:Par.nElecGps],                            ...
                 'Display',                         false,                                       ...
                 'Batch',                           false,                                       ...
                 'overwrite',                       false                                        ...
);
[Electrodes,Display,Batch,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------
nq =struct([]);

if exist(fullfile(Session.spath,[Session.name '.NeuronQuality.mat']),'file') & ~overwrite
    fullfile(Session.spath,[Session.filebase '.NeuronQuality.mat'])
    return;
end


SampleRate = Par.SampleRate;
SpkSamples = Par.SpkGrps(1).nSamples;
GoodElectrodes=[];

for el=Electrodes
    if ~exist(fullfile(Session.spath,[Session.name '.clu.' num2str(el)]),'file') continue; end
    Clu = LoadClu(fullfile(Session.spath,[Session.name '.clu.' num2str(el)]));
    if max(Clu) > 1
        GoodElectrodes = [GoodElectrodes, el];
    end
    if exist(fullfile(Session.spath,[Session.name '.spk.' num2str(el)]),'file')
        nspk = FileLength(fullfile(Session.spath,[Session.name '.spk.' num2str(el)]))/2/length(Par.ElecGp{el})/SpkSamples;
        if nspk ~= length(Clu)
            fprintf('%s - Electrode %d - length of spk file does not correspond to clu length\n\n',Session.name,el);
            nspk=-1;
        end
    else
        fprintf('%s - Electrode %d - length of spk file does not exist \n',Session.name,el);
        nspk=-1;
    end
end
%GoodElectrodes
% for non-noise cluster do the rest
%nq = struct;
for el=GoodElectrodes
    Clu = LoadClu(fullfile(Session.spath,[Session.name '.clu.' num2str(el)]));
    Res = load   (fullfile(Session.spath,[Session.name '.res.' num2str(el)]));
    
    uClu = setdiff(unique(Clu),[0 1]);
    nClu = length(uClu);
    
    if nspk>0
        %noise - comes from first 10000 spikes
        %     nNoise = sum(Clu(1:10000)==1);
        %     SpkNoise = LoadSpk([FileBase '.spk.' num2str(el)],length(Par.ElecGp{el}),SpkSamples,10000);
        %     stdSpkNoise =sq(std(SpkNoise(:,:,Clu(1:10000)==1),0,3));
        
        SampleSize = 1000;
        myClu=find(Clu==1);
        
        if length(myClu)>0
            SampleSize = min(length(myClu),SampleSize);
            RndSample = sort(myClu(randsample(length(myClu),SampleSize)));
            mySpk = LoadSpkPartial(fullfile(Session.spath,[Session.name '.spk.' num2str(el)]),...
                                   length(Par.ElecGp{el}),SpkSamples,RndSample);
            stdSpkNoise  = sq(std(mySpk,0,3));% may not need it
        else
            stdSpkNoise = 1;
        end
        %load only enough spike to sample all clusters
        % create represantative spikes sample for good cells
    end
    avSpk =[]; stdSpk = [];SpatLocal=[];SpkWidthC=[];SpkWidthL=[];SpkWidthR=[];posSpk=[];FirRate = [];AvSpkAll=[];
    leftmax=[]; rightmax=[];troughamp=[];troughSD=[];spkMaxchanAll=[];
    for cnum=1:nClu
        % get spike wavesdhapes and compute SNR
        if nspk>0 % if there was a .spk file 
            SampleSize = 1000;
            myClu=find(Clu==uClu(cnum));
            
            if length(myClu)>0
                SampleSize = min(length(myClu),SampleSize);
                RndSample = sort(myClu(randsample(length(myClu),SampleSize)));
                mySpk = LoadSpkPartial(fullfile(Session.spath,[Session.name '.spk.' num2str(el)]),...
                                       length(Par.ElecGp{el}),SpkSamples,RndSample);
                avSpk(:,:,cnum) = sq(mean(mySpk, 3)) - repmat(sq(mean(mean(mySpk, 3),2)),1,size(mySpk,2)) ;
                %remove the baseline
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%updated to here                
                stdSpk(:,:,cnum)  = sq(std(mySpk,0,3));% may not need it
                
                %find the channel of largest amp (positive or negative)
                [amps,ampch] = max(abs(avSpk(:,:,cnum)),[],2);
                [~,maxampch] = max(sq(amps));
                nch = length(Par.ElecGp{el});
                nonmax = setdiff([1:nch],maxampch);
                %compute spatial localization as ratio of max ch amplitude to the mean over all others.
                if nch>1
                    SpatLocal(cnum) = amps(maxampch)/(mean(amps(nonmax))+eps);
                else
                    SpatLocal(cnum) = 0;
                end
                myAvSpk = sq(avSpk(maxampch,:,cnum)); % largest channel spike wave for that cluster
                
                %now we need to take care of the positive spikes (we reverse them)
                minamp  = abs(min(myAvSpk));
                maxamp  = max(myAvSpk);
                %        if (minamp-maxamp)/minamp < 0.05 %(spike is more positive then negative)
                if maxamp>1.2*minamp %(spike is more positive then        negative)
                    myAvSpk = -myAvSpk; %reverse the spike
                    posSpk(cnum) = 1;
                else
                    posSpk(cnum) = 0;
                end
                
                %now let's upsample the spike waveform
                ResCoef = 10; %                                                                 WAS CHANGED ON Aug.30 2005!!!!!!!!!!
                Sample2Msec = 1000./SampleRate/ResCoef; %to get fromnew samplerate to the msec
                myAvSpk = resample(myAvSpk,ResCoef,1);
                %keyboard
                %amphalf = mean(myAvSpk)-0.5*abs(min(myAvSpk-mean(myAvSpk)));
                
                [troughamp(cnum) troughTime] = min(myAvSpk);
                pts= myAvSpk(troughTime+[-5 0 5]);
                pts=pts(:);
                troughSD(cnum) = pts'*[1 -2 1]';
                amphalf = 0.5*min(myAvSpk);
                both=0;cnt=0;halfAmpTimes=[0 0];
                while both<2
                    if myAvSpk(troughTime-cnt)>amphalf & halfAmpTimes(1)==0
                        halfAmpTimes(1)=troughTime-cnt;
                        both=both+1;
                    end
                    if myAvSpk(troughTime+cnt)>amphalf & halfAmpTimes(2)==0
                        halfAmpTimes(2)=troughTime+cnt;
                        both=both+1;
                    end
                    cnt=cnt+1;
                end
                %width
                
                %  		halfAmpTimes = LocalMinima(abs(myAvSpk-amphalf));%,5,0.5*amphalf);
                %  		halfAmpTimes = halfAmpTimes(find(halfAmpTimes-troughTime)<0.
                %  		if length(halfAmpTimes>2)
                %  			val = myAvSpk(halfAmpTimes);
                %  			[dummy ind] = sort(val);
                %  			halfAmpTimes = halfAmpTimes(ind(1:2));
                %  		end
                %  		SpkWidthC(cnum) = abs(diff(halfAmpTimes))*Sample2Msec; % this is width on the hald amplitude
                %
                SpkWidthC(cnum) = diff(halfAmpTimes)*Sample2Msec;
                
                
                %dmyAvSpk = diff(myAvSpk);
                SpkPieceR = myAvSpk(troughTime:end);
                [rightmax(cnum) SpkWidthR(cnum)] = max(SpkPieceR); % this is the distance from the trough to the rise peak
                SpkWidthR(cnum)= SpkWidthR(cnum)*Sample2Msec;
                
                SpkPieceL = myAvSpk(1:troughTime);
                SpkPieceL = flipud(SpkPieceL(:)); % to look at the right time lag
                [leftmax(cnum) SpkWidthL(cnum)] = max(SpkPieceL); % this is the distance from the peak to the trough
                SpkWidthL(cnum)= SpkWidthL(cnum)*Sample2Msec;
                
                troughTime = troughTime*Sample2Msec;
                if posSpk(cnum);	myAvSpk = -myAvSpk; end
                AvSpkAll(cnum,:) = myAvSpk;
                spkMaxchanAll(cnum,:) = maxampch;
                
                if Display
                    figure(765)
                    if cnum==2
                        clf;
                    end
                    subplotfit(cnum-1,length(uClu));
                    shift = troughTime;%SpkSamples/2*Sample2Msec;
                    plot([1:SpkSamples*ResCoef]*Sample2Msec-shift,myAvSpk);
                    axis tight
                    hold on
                    Lines(0,[],'g');%trough line
                    Lines(halfAmpTimes*Sample2Msec-shift,amphalf,'r');
                    Lines(SpkWidthR(cnum),[],'r');
                    Lines(-SpkWidthL(cnum),[],'r');
                    %Lines([-SpkWidthL SpkWidthR], troughamp,'r');
                    mystr = sprintf('El=%d Clu=%d',el,cnum);
                    title(mystr);
                end
                %keyboard
            end
            % firing rate
            myRes = Res(find(Clu==uClu(cnum)));
            if length(myRes)<3
                FirRate(cnum)=0;
            else
                ISI = diff(myRes);
                MeanISI = mean(bootstrp(100,'mean',ISI));
                FirRate(cnum) = SampleRate./MeanISI;
            end
        end
        
    end
    
    if Display
        if ~Batch
            pause
        else
            reportfig(gcf,'NeuronQuality',0,[Session.name ',El=' num2str(el)],70);
        end
        figure(765); clf
    end
    
    snr = sq(mean(mean(abs(avSpk),1),2))./mean(stdSpkNoise(:));
    %		Out = [CluNo, eDist, bRat,Refrac]

    %keyboard    
    clear mySpk SpkNoise Res Clu;
    try
        Fp = fopen(fullfile(Session.spath,[Session.name '.fet.' num2str(el)]), 'r');
        nFet = fscanf(Fp, '%d', 1);
        fclose(Fp);
        
        FQOut= FileQuality(fullfile(Session.spath,Session.name), el, SampleRate/1000*20, [1:nFet-1], SampleRate*2/1000,0);
        nq(el).eDist = FQOut(:,2);
        nq(el).bRat = FQOut(:,3);
        nq(el).Refrac = FQOut(:,4);
        %        nq(el).Clus = FQOut(:,1);
    catch
        warning('FileQuality crashed, check why');
        % uClu = unique(Clu);
        %       nq(el).Clus = uClu(:);
        nq(el).eDist = zeros(length(uClu),1);
        nq(el).bRat =  zeros(length(uClu),1);
        nq(el).Refrac =  zeros(length(uClu),1);
        
    end
    nq(el).SNR = snr;
    nq(el).SpkWidthC = SpkWidthC';
    nq(el).SpkWidthR = SpkWidthR';
    nq(el).SpkWidthL = SpkWidthL';
    %fix for empty clusters
    SpkWidthL(SpkWidthL==0)=1000000;
    nq(el).TimeSym = SpkWidthR'./SpkWidthL';
    nq(el).ElNum = repmat(el,nClu,1);
    nq(el).Clus = uClu(:);
    nq(el).IsPositive = posSpk';
    nq(el).SpatLocal=SpatLocal';
    nq(el).FirRate = FirRate';
    nq(el).AvSpk = AvSpkAll;
    nq(el).maxAmpChan = spkMaxchanAll;
    nq(el).RightMax= rightmax';
    nq(el).LeftMax= leftmax';
    nq(el).CenterMax= troughamp';
    %fix for empty clusters
    rightmax(rightmax==0)=1e6;
    leftmax(leftmax==0)=1e6;
    nq(el).AmpSym = (abs(rightmax)'-abs(leftmax)')./(abs(rightmax)'+abs(leftmax)');
    nq(el).troughSD = troughSD';
    %keyboard
    
end


nq =CatStruct(nq);
save(fullfile(Session.spath,[Session.name '.NeuronQuality.mat']),'nq');


% END MAIN -----------------------------------------------------------------------------------------

