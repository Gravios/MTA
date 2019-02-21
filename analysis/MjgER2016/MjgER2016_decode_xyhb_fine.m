% MjgER2016_decode_xyhb_fine
%
% Bayesian decoding of position and behavioral features
% 
% 
%



MjgER2016_load_data();

if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
numComp = size(eigVec{1},2);
pfindex = 1;
MjgER2016_load_bhv_erpPCA_scores();
% output:
%    fsrcz
%    FSrC
%    rmaps
clear('pfd','tags','eigVec','eigVar','eigScore','validDims','zrmMean','zrmStd',...
      'clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','rmapsShuffledMean',...
      'rmapsShuffled','FSrM','FSrS','fsrsMean','fsrsStd','fsrsMean','fsrsStd');


for trialIndex = 17:23;    
Trial = Trials{trialIndex}; 
unitSubset = units{trialIndex};
[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyh',[],[],0.03,true);
end


trialIndex = 20;    
Trial = Trials{trialIndex}; 
unitSubset = units{trialIndex};

sampleRate = 250;   % Hz
spikeWindow = 0.025; % ms
mode = 'xyhi'; % alt vals: 'xy'

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit','ripple'};

% LOAD state collection
% LOAD subject position object
stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
fxyz = filter(copy(xyz),'ButFilter',3,20,'low');

% LOAD lfp
% COMPUTE theta LFP phase
lfp = load(Trial,'lfp',72);
%lfp = load(Trial,'lfp',84);
phz = lfp.phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);    

% COMPUTE speed
vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

% LOAD pitch and 
fet = fet_HB_pitchB(Trial,sampleRate);
tfet = copy(fet);
fet.data = fet(:,1);
ffet = filter(copy(fet),'ButFilter',3,1,'low');

tvec = circshift(fxyz(:,'hcom',[1,2]),-round(sampleRate.*0.05))-fxyz(:,'hcom',[1,2]);
tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);
hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
bvec = fxyz(:,'spine_lower',[1,2])-fxyz(:,'spine_upper',[1,2]);
bvec = sq(bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3))));
bvec = cat(3,bvec,sq(bvec)*[0,-1;1,0]);

% CONVERT MTAStateCollection into a state matrix 
stcm = stc2mat(stc,xyz,states);


ufr = load(Trial,'ufr',xyz,[],unitSubset,spikeWindow,true,'gauss');    
unitInclusion = sum(ufr.data>0.2,2);

[posEstCom,posEstMax,posteriorMax] = bhv_decode(Trial,sampleRate,unitSubset,mode,[],[],spikeWindow,true);

[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyh',[],[],0.025,false);
[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyhi',[],[],0.025,false);
[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyb',[],[],0.025,false);
[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyz',[],[],0.025,false);


[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xyhi',[],[],0.125,false);

ind = unitInclusion>=4 & any(stcm==1,2) & posteriorMax>0.001;
sum(ind)

%ind = ':';

%derror = posEstMax;
derror = posEstSax;
%derror = posEstCom;
figure()
subplot(511);  plot(derror(ind,[1]));  %Lines([0;cumsum(diff([rper.data,1,2))],[],'k');
hold('on');    plot(xyz(ind,'hcom',1),'k');
subplot(512);  plot(derror(ind,[2]));  %Lines([0;cumsum(diff(rper.data,1,2))],[],'k');
hold('on');    plot(xyz(ind,'hcom',2),'k');
subplot(513);  plot(derror(ind,[3]));  %Lines([0;cumsum(diff([rper.data,1,2))],[],'k');
hold('on');    plot(fet(ind,1),'k');
%subplot(514);  plot(derror(ind,[4]));  %Lines([0;cumsum(diff(rper.data,1,2))],[],'k');
%hold('on');    plot(fet(ind,2),'k');
subplot(515);  %plot(mufr(ind));  %Lines([0;cumsum(diff(rper.data,1,2))],[],'k');
hold('on');    plot(posteriorMax(ind)*10);
linkaxes(findobj(gcf,'Type','Axes'),'x');
clear('derror');


%% some plots of stuff
ind = ':';
decError = [multiprod(posEstSax(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]),fet(ind,1)-posEstSax(ind,3)];
decError(~(unitInclusion>=2 & any(stcm==1,2) & posteriorMax>0.0005)) = nan;

figure,
hold('on');
plot(ang(:,'head_back','head_front',1));
plot(atan2(decError(:,2),decError(:,1)),'.r');
plot(phz(:,1),'c')

% Decompose decoding error by theta phase and state in egocentric frame of reference
figure();
ns = 6;
decError = [multiprod(posEstCom(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]),fet(:,1)-posEstCom(:,3)];
%decError = [multiprod(posEstSax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]),fet(:,1)-posEstSax(:,3)];
for s = 1:ns; %states
    cind = unitInclusion>=4 & posteriorMax>0.001;
    sind =      stcm(:,1)   ...
        &  any(stcm(:,s),2);

    ind = cind&sind;

    % Forward decoding projection onto head
    subplot2(ns,3,s,1);
    out = hist2([[decError(ind,1);decError(ind,1)],[phz(ind,1);phz(ind,1)+2*pi]],...
          linspace(-300,300,50),linspace(-pi,pi*3,50));
    imagesc(linspace(-300,300,50),linspace(-pi,pi*3,50),bsxfun(@rdivide,out,max(out))');
    axis('xy');
    title('forward projection')
    xlim([-300,300]);
    xlabel({states{s},'mm'});
    ylabel('theta phase');
    % Lateral decoding projection onto head
    subplot2(ns,3,s,2);
    out = hist2([[decError(ind,2);decError(ind,2)],[phz(ind,1);phz(ind,1)+2*pi]],...
          linspace(-300,300,50),linspace(-pi,pi*3,50));
    imagesc(linspace(-300,300,50),linspace(-pi,pi*3,50),bsxfun(@rdivide,out,max(out))');
    axis('xy');    
    title('lateral projection')
    xlabel('mm');
    ylabel('theta phase');

    % Pitch decoding error
    subplot2(ns,3,s,3);
    out = hist2([[decError(ind,3);decError(ind,3)],[phz(ind,1);phz(ind,1)+2*pi]],...
          linspace(-pi/3,pi/3,75),linspace(-pi,pi*3,50));
    imagesc(linspace(-pi/3,pi/3,75),linspace(-pi,pi*3,50),bsxfun(@rdivide,out,max(out))');
    axis('xy');    
    title('head pitch projection');
    xlabel('rad');
    ylabel('theta phase');
end


figure
%decError = [multiprod(posEstCom(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]),fet(:,1)-posEstCom(:,3)];
decError = [multiprod(posEstSax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]),fet(:,1)-posEstSax(:,3)];
%decError = [multiprod(posEstMax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]),fet(:,1)-posEstMax(:,3)];
cind = unitInclusion>=4 & posteriorMax>0.001;
sind =      stcm(:,1)   ...
    &  (stcm(:,3)==3|stcm(:,5)==5);
ind = cind&sind;
ind = ind&circ_dist(circshift(ffet(:,1),-1),ffet(:,1))>0.0001;
subplot(121);
hist2([[phz(ind,1);phz(ind,1)+2*pi],[decError(ind,3);decError(ind,3)]],linspace(-pi,pi*3,50),linspace(-pi/3,pi/3,75));
ind = cind&sind;
ind = ind&circ_dist(circshift(ffet(:,1),-1),ffet(:,1))<-0.0001;
subplot(122);
hist2([[phz(ind,1);phz(ind,1)+2*pi],[decError(ind,3);decError(ind,3)]],linspace(-pi,pi*3,50),linspace(-pi/3,pi/3,75));


figure,plot([fet.data,posEstSax(:,3).*double(unitInclusion>=3)])


s = 1:6;
cind = unitInclusion>=3 & posteriorMax>0.0001;
sind =      stcm(:,1)   ...
         &  stcm(:,s)   ...    
         & ~stcm(:,2)   ...
         & ~stcm(:,7)   ...          
         & ~stcm(:,8);

ind = cind&sind;
sum(ind)


% plot the 

hfig = figure();
% STATE - locomotion
nPart = 6
speedInd = 2

clf();
hfig.Units = 'centimeters';
hfig.Position = [0,0,30,5*nPart];
for v = 1:nPart
    y = nPart-v+1;
    % select data with in speed partition v
    ind = cind  &  sind  ...
                &  vxy(:,speedInd)>modelTraj.partitions(v)   ...
                &  vxy(:,speedInd)<modelTraj.partitions(v+1);
% REFERENCE trajectory coordinate system
    %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
    decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
    
    % Anterior-Posterior axis
    subplot2(nPart,6,y,1);
    % JPDF phase X error
    hist2([[decError(:,1);decError(:,1)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])-2*pi,'Color','m');    
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300]),     'Color','m');
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])+2*pi,'Color','m');
    title(['rho: ',num2str(modelTraj.rho(v))]);
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% medial-lateral axis
    subplot2(nPart,6,y,2);
    hist2([[decError(:,2);decError(:,2)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% direction
    subplot2(nPart,6,y,3);
    hist2([[atan2([decError(:,2);decError(:,2)],...
                  [decError(:,1);decError(:,1)])],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-pi,pi,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('yaw');
    ylabel('theta phase');
    Lines(0,[],'k');
% REFERENCE head coordinate system
    %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
    decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
    subplot2(nPart,6,y,4);
    hist2([[decError(:,1);decError(:,1)],...
          [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])-2*pi,'Color','m');
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300]),'Color','m');
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])+2*pi,'Color','m');
    
    title(['rho: ',num2str(modelHead.rho(v))]);
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% lateral
    subplot2(nPart,6,y,5);
    hist2([[decError(:,2);decError(:,2)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% yaw
    subplot2(nPart,6,y,6);
    hist2([atan2([decError(:,2);decError(:,2)],...
                 [decError(:,1);decError(:,1)]),...
          [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-pi,pi,40),...
          linspace(-pi,pi*3,30));
    Lines(0,[],'k');
    xlabel('yaw');
    ylabel('theta phase');
end



ForAllSubplots('hax=gca;hax.Units=''centimeters'';hax.Position(3:4)=[3,2];')

af(@(hax) set(hax,'Position',hax.Position+[0,-2,0,0]), findobj(gcf,'Type','Axes'));

axPos = cell2mat([get(findobj(gcf,'Type','Axes'),'Position')]);
axYPos = unique(axPos(:,2));

fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);
Lines(hfig.Position(3)/2+0.2,[],'k');

Lines([],axYPos+3,'k');

for v = 1:nPart, 
    text(1,axYPos(nPart-v+1)+1,...
         {[tag,' speed:'],...
          num2str(modelTraj.partitions(nPart-v+2)),...
          '',...
          num2str(modelTraj.partitions(nPart-v+1))},...
         'Rotation',0);
end

ht = text(hfig.Position(3).*0.3,hfig.Position(4)*.9,...
          'Head Movement Basis', 'HorizontalAlignment','center');
ht = text(hfig.Position(3).*0.7,hfig.Position(4)*.9,...
          'Head Direction Basis','HorizontalAlignment','center');

text(1,hfig.Position(4)*.95,...
     {[Trial.filebase,': locomation and pause'],...
      ['unitInclusion >= 2'],...
      [ 'posteriorMax >  0.003'],...
      ['JPDF of projected decoded positions and theta phase,'],...
      ['partition over equal partitions of ' tag ' speed']})

print(hfig,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_decodedXthetaPhaseX',tag,'Speed_',mode,'_',Trial.filebase,'.eps']]);
print(hfig,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
      ['dc_decodedXthetaPhaseX',tag,'Speed_',mode,'_',Trial.filebase,'.png']]);

% SAVE vars
    save(fullfile(Trial.spath,...
                  [Trial.filebase,'_dc_decodedXthetaPhaseX',tag,'Speed_',mode,'.mat']),...
         'modelTraj',...
         'modelHead');








% TEMPORAL shift analysis
% find point in trajectory closest to decoded point
% report time shift and distance

shiftvals = -round(sampleRate*2):round(sampleRate/50):round(sampleRate*2);
% SET time bins
slidingTimeWindowSize = round(sampleRate*1);
timevals = 1:size(xyz,1);
timevals(~nniz(xyz)) = [];
timevals = timevals(1:slidingTimeWindowSize:numel(timevals));
timeBinInds  = discretize([1:size(xyz,1)]',timevals);
badTimeBinInds = find(diff(timevals)>slidingTimeWindowSize)+1;
timeBinInds(ismember(timeBinInds,badTimeBinInds)) = 0;
% SET phase bins
phasevals = linspace(-pi,pi,13);
phaseBinInds = discretize(phz.data,phasevals);
% SET speed bins
speedvals = linspace(0,1.8,7);
speedBinInds = discretize(vxy(:,2),speedvals);

eposTPRErrorXY = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);
%eposTPRErrorHP = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);
%eposTPRErrorBP = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);

eposTPRErrorXYTrj = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1,2]);

meanTPRSpeed = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPostMax = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPhase = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);

for time = 1:numel(timevals)-2,
    %tic
    disp([num2str(time),' of ',num2str(numel(timevals)-3)])

    timeWindow = time-2:time+2;
    
    timeWindowIndex = ismember(timeBinInds,timeWindow);

    txyz = sq(xyz(timeWindowIndex,'nose',[1,2]));
    tvxy = vxy(timeWindowIndex,2);
    tbvec = bvec(timeWindowIndex,:,:);    
    tposteriorMax = posteriorMax(timeWindowIndex);
    tphz = phz(timeWindowIndex);    
    
    tind =   unitInclusion(timeWindowIndex)>=3  ...
             & stcm(timeWindowIndex,1)==1       ...
             & stcm(timeWindowIndex,2)~=2       ...
             & timeBinInds((timeWindowIndex))==time...
             & tposteriorMax>0.002;

    tempPhaseBinInds = phaseBinInds(timeWindowIndex);
    tempSpeedBinInds = speedBinInds(timeWindowIndex);

    eposThetaPhaseRes = nan([size(tind,1),2]);
    eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstCom(timeBinInds==time,1:2);
    %eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstMax(timeBinInds==time,1:2);
% $$$     eposThetaPhaseRes = nan([size(tind,1),size(posEstCom,2)]);
% $$$     eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstCom(timeBinInds==time,1:2);

    for speed = 1:numel(speedvals)-1,   
        for phase = 1:numel(phasevals)-1,   
            ind = tind  &  tempPhaseBinInds==phase  &  tempSpeedBinInds==speed;
            
            if sum(ind)>3
                meanTPRSpeed(time,phase,speed) = mean(tvxy(ind));
                meanTPRPostMax(time,phase,speed) = mean(tposteriorMax(ind));
                meanTPRPhase(time,phase,speed) = circ_mean(tphz(ind));
                tepos = eposThetaPhaseRes;
                tepos(~ind,:) = nan;


                                    
                for shift = 1:numel(shiftvals)
                    eposTPRErrorXYTrj(time,shift,phase,speed,:) = ...
                       mean(multiprod(circshift(tepos(:,1:2),shiftvals(shift))-txyz,tbvec,2,[2,3]),'omitnan');
                    
                    eposTPRErrorXY(time,shift,phase,speed) = ...
                        mean(sqrt(sum((txyz-circshift(tepos,shiftvals(shift))).^2,2)),'omitnan');

% $$$                     %                    eposTPRErrorHP(time,shift,phase,speed) = mean(sqrt(sum((fet(:,1) ...
% $$$                                        -circshift(eposThetaPhaseRes(:,3),shiftvals(shift))).^2,2)),'omitnan');
% $$$                     %eposTPRErrorBP(time,shift,phase,speed) = mean(sqrt(sum((fet(:,2)...
% $$$                                        -circshift(eposThetaPhaseRes(:,4),shiftvals(shift))).^2,2)),'omitnan');
                end
            end
        end
    end
    %toc
end



% $$$ figure();
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ subplot2(2,numel(speedvals)-1,1,speed);    
% $$$     imagesc(phasevals,shiftvals/sampleRate,eposTPRErrorXY(:,:,speed));axis('xy');
% $$$ caxis([100,300])
% $$$ subplot2(2,numel(speedvals)-1,2,speed);
% $$$     imagesc(phasevals,shiftvals/sampleRate0,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min( ...
% $$$         eposTPRErrorXY(:,:,speed)))+eps));axis('xy');
% $$$     title(num2str(mean(10.^speedvals(speed:speed+1))))
% $$$ end



derror = eposTPRErrorXY;
%derror = abs(eposTPRErrorXYTrj(:,1));
% $$$ derror = eposTPRErrorBP;
% $$$ derror = eposTPRErrorHP;

mind = [];
minv = [];
%for time = 1:numel(timevals)-1,
for t = 1:time
    for speed = 1:numel(speedvals)-1,
        %[~,mind(t,speed,:)] = max(RectFilter(RectFilter(1./(bsxfun(@rdivide,sq(derror(t,:,:,speed)),min(sq(derror(t,:,:,speed))))+eps)',3,1)',3,1));
        [minv(t,speed,:),mind(t,speed,:)] = min(RectFilter(RectFilter(sq(derror(t,:,:,speed))',3,1)',3,1));
    end
end
mind(mind==0) = nan;
minv(minv==0) = nan;


shiftTimeBins = shiftvals(2:end)./sampleRate;

figure;
for s = 1:numel(speedvals)-1,   
    for p = 1:numel(phasevals)-1,   
        subplot2(numel(phasevals)-1,numel(speedvals)-1,p,s);

        ind = mind(:,s,p)~=1&mind(:,s,p)~=201;

        bar(shiftTimeBins,histc(shiftTimeBins(mind(ind,s,p)),shiftTimeBins),'histc');
        %hist2([shiftTimeBins(mind(ind,s,p))',minv(ind,s,p)],linspace(-1,2,30),linspace(40,800,6));
% $$$         hist2([shiftTimeBins(mind(ind,s,p))',minv(ind,s,p)],linspace(-1,2,30),linspace(40,800,6));
        Lines(median(shiftTimeBins(mind(ind,s,p)),'omitnan'),[],'g');
        Lines(0,[],'m');

% $$$         scatter(log10(eposTPRErrorXY([ind;false;false],50,p,s)),log10(minv(ind,s,p)),5, ...
% $$$                 shiftTimeBins(mind(ind,s,p)),'filled');
% $$$         xlim([1,2.7]);
% $$$         ylim([1,2.7]);
        
% $$$         scatter(log10([1;1;eposTPRErrorXY([ind;false;false],50,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         xlim([1,2.7]);
% $$$         ylim([1,2.7]);
% $$$         scatter(([0;0;eposTPRErrorXYTrj([ind;false;false],50,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         xlim([-400,400]);
% $$$         ylim([1,2.7]);
        
        xlim(shiftvals([1,end])./sampleRate);
        ylim([0,30]);
        if s==1,ylabel(num2str(mean(phasevals([p:p+1]))));end
        if p==(numel(phasevals)-1),xlabel(num2str(speedvals([s:s+1])));end
            
    end
end




% $$$ 
% $$$ cm = reshape(permute(repmat(jet(numel(shiftvals)-1),[1,1,20]),[1,3,2]),[],3);
% $$$ figure();
% $$$ scatter(reshape(repmat(phasevals(1:end-1),[numel(shiftvals)-1,1]),[],1),shiftvals(reshape(mind,[],1))/40,20,cm,'filled');

figure();
imagesc(phasevals,speedvals,(mind-40)/40)
axis('xy');
caxis([-0.4,0.4]);

imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,18),min(eposTPRErrorXY(:,:,18)))+eps));axis('xy');

figure();
subplot(131);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY,min(eposTPRErrorXY))+eps));axis('xy');
subplot(132);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorHP,min(eposTPRErrorHP))+eps));axis('xy');
subplot(133);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorBP,min(eposTPRErrorBP))+eps));axis('xy');







% OLD 
% $$$ 
% $$$ figure();
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ subplot2(2,numel(speedvals)-1,1,speed);    
% $$$     imagesc(phasevals,shiftvals/40,eposTPRErrorXY(:,:,speed));axis('xy');
% $$$ caxis([100,300])
% $$$ subplot2(2,numel(speedvals)-1,2,speed);
% $$$     imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min( ...
% $$$         eposTPRErrorXY(:,:,speed)))+eps));axis('xy');
% $$$     title(num2str(mean(10.^speedvals(speed:speed+1))))
% $$$ end
% $$$ 
% $$$ 
% $$$ mind = [];
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ [~,mind(speed,:)] = max(RectFilter(RectFilter(1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min(eposTPRErrorXY(:,:,speed)))+eps)',3,1)',3,1));
% $$$ end
% $$$ % $$$ 
% $$$ % $$$ cm = reshape(permute(repmat(jet(numel(shiftvals)-1),[1,1,20]),[1,3,2]),[],3);
% $$$ % $$$ figure();
% $$$ % $$$ scatter(reshape(repmat(phasevals(1:end-1),[numel(shiftvals)-1,1]),[],1),shiftvals(reshape(mind,[],1))/40,20,cm,'filled');
% $$$ 
% $$$ figure();
% $$$ imagesc(phasevals,speedvals,(mind-40)/40)
% $$$ axis('xy');
% $$$ caxis([-0.4,0.4]);
% $$$ 
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,18),min(eposTPRErrorXY(:,:,18)))+eps));axis('xy');
% $$$ 
% $$$ figure();
% $$$ subplot(131);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY,min(eposTPRErrorXY))+eps));axis('xy');
% $$$ subplot(132);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorHP,min(eposTPRErrorHP))+eps));axis('xy');
% $$$ subplot(133);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorBP,min(eposTPRErrorBP))+eps));axis('xy');






    


% CREATE video of decoding
% $$$ xbinsReal = discretize(xyz(:,'nose',1),binsE{1});
% $$$ ybinsReal = discretize(xyz(:,'nose',2),binsE{2});
% $$$ 
% $$$ xbinsReal = discretize(tpos(:,1),binsE{1});
% $$$ ybinsReal = discretize(tpos(:,2),binsE{2});
% $$$ hbinsReal = discretize(tpos(:,3),binsE{3});
% $$$ bbinsReal = discretize(tpos(:,4),binsE{4});
% $$$ 
% $$$ 
% $$$ [C,H] = bhv_contours();
% $$$ 
% $$$ H{1}.ZData(sqrt(sum(cat(3,H{1}.XData,H{1}.YData).^2,3))<0.2)=0;
% $$$ 
% $$$ 
% $$$ vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/decode/posterior_example_immobile.avi','Uncompressed AVI');
% $$$ open(vidObj);
% $$$ 
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ pause(0.1);
% $$$ hfig.Position(3:4) = [30,20];
% $$$ pause(0.1);
% $$$ axPosition = axes('Units','centimeters','Position',[2,4,6,6]);
% $$$ hold(axPosition,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axPosture = axes('Units','centimeters','Position',[2,12,6,6]);
% $$$ hold(axPosture,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axBg = axes('Position',[0 0 1 1],'Visible','off');    
% $$$ 
% $$$ axPosX = axes('Units',         'centimeters',...
% $$$               'Position',      [13,16,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('X (cm)');
% $$$ haxLinesTime(1) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,13,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('Y (cm)');
% $$$ haxLinesTime(2) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,10,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Head','Pitch (rad)'});
% $$$ haxLinesTime(3) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,7,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Body','Pitch (rad)'});
% $$$ haxLinesTime(4) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axStc  = axes('Units',         'centimeters',...
% $$$               'Position',      [13, 4,14,3],...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on');
% $$$ plotSTC(stc,1,[],fliplr({'rear','loc','pause'}),fliplr('rbk'));
% $$$ axStc.YTickLabel = {'pause','loc','rear'};
% $$$ xlim([310,375]);
% $$$ xlabel('Time (s)');
% $$$ haxLinesTime(5) = animatedline([310,310],[1,4],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ limitsBound = [-50,50;-50,50;-1.5,1.5;-1.5,1.5,;1,4];
% $$$ 
% $$$ for t = 3100:1:3750
% $$$ 
% $$$     axes(axPosition);    
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-50,50,100),linspace(-50,50,100),zeros([100,100]));
% $$$         axPositionCirc = circle(0,0,42,'--r');        
% $$$     end
% $$$ 
% $$$     if ~isnan(hbinsReal(t))&~isnan(bbinsReal(t)),
% $$$         axPositionIm = imagesc(binsE{1}/decodingSampleRate,binsE{2}/decodingSampleRate,sq(tE(:,:,hbinsReal(t),bbinsReal(t),t))');
% $$$         axPositionScXyz = scatter(xyz(t,'nose',1)/decodingSampleRate,xyz(t,'nose',2)/decodingSampleRate,40,'g','filled');
% $$$         axPositionScTpos = scatter(tpos(t,1)/decodingSampleRate,tpos(t,2)/decodingSampleRate,40,'m','filled');
% $$$     end    
% $$$     uistack(axPositionCirc,'top')
% $$$ 
% $$$     if t == 3100,        
% $$$ 
% $$$     end
% $$$     
% $$$     xlabel('X Position (cm)');
% $$$     ylabel('Y Position (cm)');
% $$$     title('Decoded Position');
% $$$     xlim([-50,50]);    
% $$$     ylim([-50,50]);    
% $$$     daspect([1,1,1]);
% $$$ 
% $$$     axes(axPosture);
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-1.65,1.4,100),linspace(-1,1.75,100),zeros([100,100]));
% $$$     end
% $$$     
% $$$     if ~isnan(xbinsReal(t))&~isnan(ybinsReal(t)),
% $$$         axPostureIm = imagesc(binsE{3},binsE{4},sq(tE(xbinsReal(t),ybinsReal(t),:,:,t))');
% $$$         axPostureScFet = scatter(fet(t,1),fet(t,2),40,'g','filled');
% $$$         axPostureScTpos = scatter(tpos(t,3),tpos(t,4),40,'m','filled');
% $$$     end
% $$$     xlabel('Head-Body Pitch (rad)');
% $$$     ylabel('Body Pitch (rad)')
% $$$     title('Decoded Posture');
% $$$     xlim([-1.65,1.4]);
% $$$     ylim([-1,1.75]);
% $$$     
% $$$     if t == 3100,
% $$$         HC = cf(@(h) copyobj(h,gca), H);
% $$$         h = plot(0,0,'w');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'c');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'r');
% $$$         h.Visible = 'off';        
% $$$         lax = legend({'real','decode','low','high','rear'},'Location','eastoutside');
% $$$         lax.Position = lax.Position+[0.1,0,0,0];
% $$$         lax.Color = [0.6,0.6,0.6];
% $$$     end
% $$$     cf(@(h) uistack(h,'top'), HC);
% $$$     
% $$$ 
% $$$     
% $$$     axes(axBg);
% $$$     cla(axBg);
% $$$     text(0.14,0.1,['time: ',num2str(t/decodingSampleRate),' s']);
% $$$     
% $$$     af(@(h) clearpoints(h), haxLinesTime);
% $$$     for h = 1:numel(haxLinesTime),
% $$$          addpoints(haxLinesTime(h),[t/decodingSampleRate,t/decodingSampleRate],limitsBound(h,:));
% $$$     end
% $$$ 
% $$$     drawnow('update');
% $$$ 
% $$$     writeVideo(vidObj,getframe());
% $$$     
% $$$     % Clear axes
% $$$     delete([axPositionIm,axPositionScXyz,axPositionScTpos]);    
% $$$     delete([axPostureIm,axPostureScFet,axPostureScTpos]);
% $$$     
% $$$ end
% $$$ close(vidObj);

% $$$ 
% $$$ % DIAGNOSTIC plot of decoded and real positions in physical and behavioral spaces
% $$$ figure();
% $$$ subplot(511);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10);
% $$$ title('x pos');
% $$$ subplot(512);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10);
% $$$ title('y pos');
% $$$ subplot(513);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1));
% $$$ title('head pitch');
% $$$ subplot(514);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2));
% $$$ title('body pitch');
% $$$ subplot(515);
% $$$ plotSTC(stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause','theta'}),fliplr('rbcgkm'));
% $$$ xlabel('Time (s)');
% $$$ linkaxes(get(gcf,'Children'),'x');
% $$$ linkaxes(get(gcf,'Children'),'x');


% ENCAPSULATE decoded position in MTAData object

rpos = MTADfet.encapsulate(Trial,                                                             ... % MTATrial object
                           tpos,                                                              ... % Data
                           decodingSampleRate,                                                ... % Sample rate
                           ['bayesiandecoded_',pfs.parameters.type,'_',pfs.parameters.states],... % Name
                           ['bd',pfs.parameters.type],                                        ... % Label
                           'b'                                                                ... % key
);

ind = [stc{'a&t-m-s',xyz.sampleRate}];


% FIGURE - jpdf of physical and behavioral spaces
figure,
subplot(121);
out = hist2(rpos(ind,[1,2])/10,linspace([-50,50,50]),linspace([-50,50,50]));
imagesc(linspace([-50,50,50]),linspace([-50,50,50]),out'./sum(out(:)));
axis('xy');
caxis([0,max(caxis)/2]);
daspect([1,1,1]);
title({'Decoded Position',Trial.filebase});
xlabel('cm');
ylabel('cm');
if strcmp(pfstype,'xyhb');
    subplot(122);
    out = hist2(rpos(ind,[3,4]),linspace(-1.5,1,50),linspace(-1,1.8,50));
    imagesc(linspace(-1.5,1,50),linspace(-1,1.8,50),out'./sum(out(:)));
    axis('xy');
    caxis([0,max(caxis)/2]);
    daspect([1,1,1]);
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.png']]);



meanPopRate = sum(ufr(ind,:),2);
meanPopIncl = sum(ufr(ind,:)>1,2);

rposError = [sqrt((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2)];
spatialError = sqrt(sum((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2,2));
if strcmp(pfstype,'xyhb'),
    rposError = cat(2,rposError,sqrt((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2));  
    bhvError = sqrt(sum((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2,2));
else
    bhvError = [];
    bhvErrorCondPos = [];
end






% FIGURE - jpdf of phisical and behavioral errors
if strcmp(pfstype,'xyhb'),
figure();
hist2([spatialError/10,bhvError],linspace(0,20,25),linspace(0,1.4,25));
title({'Spatial Vs Behavioral Error',[Trial.filebase,' ',num2str(numel(unitSubset)),' units']});
xlabel('cm');
ylabel('rad');
hax = gca();
hax.Units = 'centimeters';
hax.Position = [hax.Position([1,2]),4,4];
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.png']]);
end


% COMPUTE mean error for decoded components individually
shuffle = @(x) circshift(x,randi(size(x,1)));
ind = [stc{'a&t-m-s',xyz.sampleRate}];
meanPopRate = sum(ufr(ind,:),2);
txss = sq(xyz(ind,'nose',[1,2]));
tfss = fet(ind,[1,2]);
trps = rpos(ind,:);
meanError = mean(rposError);
% COMPUTE shuffled error for decoded components
meanErrorShuffled = nan([1000,numel(pfsBins)]);
bhvErrorShuffled = [];
mbes = [];
for i = 1:1000,
    sptErrorShuffled = sqrt((shuffle(txss)-trps(:,[1,2])).^2);
    if strcmp(pfstype,'xyhb'),
        bhvErrorShuffled = sqrt((shuffle(tfss)-trps(:,[3,4])).^2);
        mbes = mean(bhvErrorShuffled(meanPopRate>0.6,:));
    end
    meanErrorShuffled(i,:) = [mean(sptErrorShuffled(meanPopRate>0.6,:)), mbes];
end

% COMPUTE zscore of
for e =1:numel(meanError),
zsMSError(e) = (meanError(e)-mean(meanErrorShuffled(:,e)))/std(meanErrorShuffled(:,e));
end



% FIGURE - mean error conditiond on position
posBins = linspace(-50,50,25);
xBins = discretize(xyz(ind,'nose',1)./10,posBins);
yBins = discretize(xyz(ind,'nose',2)./10,posBins);
aind = xBins>0&yBins>0;
figure();
subplot(121);
spatialErrorCondPos = accumarray([xBins(aind),yBins(aind)],spatialError(aind)/10,[numel(posBins),numel(posBins)],@median);
spatialErrorCondPos(spatialErrorCondPos==0) = nan;
imagescnan({posBins,posBins,spatialErrorCondPos'},[0,20],'linear',true,'colorMap',@jet);axis('xy');
title({'Mean Decoded Spatial Error','Conditioned on Position'});
if strcmp(pfstype,'xyhb'),
    subplot(122);
    bhvErrorCondPos = accumarray([xBins(aind),yBins(aind)],bhvError(aind),[numel(posBins),numel(posBins)],@median);
    bhvErrorCondPos(bhvErrorCondPos==0) = nan;
    imagescnan({posBins,posBins,bhvErrorCondPos'},[0,1],'linear',true,'colorMap',@jet);axis('xy');
    title({'Mean Decoded Pitch Error','Conditioned on Position'});
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
hax(1).Position = [hax(1).Position([1,2]),4,4];
hax(2).Position = [hax(2).Position([1,2]),0.5,4];
if strcmp(pfstype,'xyhb'),
    hax(3).Position = [hax(3).Position([1,2]),4,4];
    hax(4).Position = [hax(4).Position([1,2]),0.5,4];
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.png']]);



% FIGURE - Independence of errors
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,14,16];

subplot(423);hist2(mud([rposError(:,[1,2])]),50,50);title('x vs y');
if strcmp(pfstype,'xyhb'),
    subplot(421);hist2(mud([spatialError,bhvError]),50,50);title('xy vs hb');
    subplot(424);hist2(mud([rposError(:,[3,4])]),50,50);title('h vs b');
    subplot(425);hist2(mud([rposError(:,[1,3])]),50,50);title('x vs h');
    subplot(426);hist2(mud([rposError(:,[2,4])]),50,50);title('y vs b');
    subplot(427);hist2(mud([rposError(:,[1,4])]),50,50);title('x vs b');
    subplot(428);hist2(mud([rposError(:,[2,3])]),50,50);title('y vs h');
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]),  hax);
axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
text(0.5,0.9,{Trial.filebase,['units: ',num2str(numel(unitSubset))]});
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.png']]);





save(fullfile(MTA_PROJECT_PATH,'analysis',['bhv_decode_',pfstype,'_',Trial.filebase,'.mat']),...
     'ind',                         ...
     'rposError',                   ...
     'spatialError',                ...
     'bhvError',                    ...
     'spatialErrorCondPos',         ...     
     'bhvErrorCondPos',             ...
     'posBins',                     ...
     'meanError',                   ...
     'meanErrorShuffled',           ...
     'zsMSError',                   ...
     'meanPopRate',                 ...
     'meanPopIncl',                 ...
     'unitSubset'                   ...
);

end


clear('ds');
trialSubset = [1:4,6,7,17,18,20:23];
pfstypes = {'xy','xyhb'};
for trialIndex =  1:numel(trialSubset);
    Trial = Trials{trialSubset(trialIndex)}; 
    for typeIndex = 1:numel(pfstypes)
        ds(trialIndex,typeIndex) = load(fullfile(MTA_PROJECT_PATH,'analysis',...
                              ['bhv_decode_',pfstypes{typeIndex},'_',Trial.filebase,'.mat']));
    end
end

me2d = reshape([ds(:,1).meanError],[2,12])';
me4d = reshape([ds(:,2).meanError],[4,12])';

me2d = [ds(:,1).spatialError']';
me4d = reshape([ds(:,2).meanError],[4,12])';



me2dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,1)));
me4dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,2)));

me2dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,1)));
me4dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,2)));

me2s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,1)));
me4s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,2)));

me4p = cell2mat(af(@(s) median(s.bhvError(:,1)), ds(:,2)));

nunits = cell2mat(af(@(s) numel(s.unitSubset(:)), ds(:,1)));

figure();
subplot
hold('on');
scatter(nunits,me2dx,20,'b','filled')
scatter(nunits,me4dx,20,'r','filled')

figure();
hold('on');
scatter(nunits,me2dy,20,'b','filled')
scatter(nunits,me4dy,20,'r','filled')



figure();
subplot(121);
hold('on');
scatter(nunits,me2s/10,20,'b','filled');
scatter(nunits,me4s/10,20,'r','filled');
xlabel('Number of Units');
ylabel('Distance (cm)');
title('Median Spatial Error');
legend({'2d place fields','4d place fields'},'location','eastoutside');
grid('on')
set(gca,'XTick',[0,25,50,75,100]);
set(gca,'YTick',[4,5,6,7,8,9,10]);

subplot(122);
scatter(nunits,me4p,20,'r','filled')
xlabel('Number of Units');
ylabel('Distance (rad)');
title('Median Behavioral Error');
legend({'4d place fields'},'location','eastoutside');
grid('on');
set(gca,'XTick',[0,25,50,75,100]);

hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),3,3]),  hax);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.png']]);




% 20180725 ----------------------------------------------
trialIndex = 20;
Trial = Trials{trialIndex};
stc = Trial.load('stc','msnn_ppsvd_raux');
unitSubset = units{trialIndex};

xyz = preproc_xyz(Trial,'trb');
xyz.resample(250);

vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','head_back'});

states = {'rear','hloc','hpause','lloc','lpause','groom','sit','theta','spw'};

drz = xyz.copy();
drz.data = compute_drz(Trial,unitSubset,pfs,'feature',xyz.copy());

lfp = Trial.load('lfp',[68,72,76,82]);

phz = phase(resample(copy(lfp),xyz),[5,12]);

spk = Trial.load('spk',xyz.sampleRate,[],unitSubset);

ufrAll = Trial.load('ufr',xyz,[],[],0.04);

unitsInt = select_units(Trial,'int');
ufrInt = Trial.load('ufr',xyz,[],unitsInt,0.04);

int = Trial.load('spk',xyz.sampleRate,'',unitsInt);

figure
for unit = 1:numel(unitsInt),
    subplotfit(unit,20);
    rose(phz(int(unitsInt(unit)),1));
end



ufrIntObj = ufrInt.copy();
ufrIntObj.data = mean(ufrInt(:,[5,6,10]),2);
uiphz = ufrIntObj.phase([5,12]);

ufrAllObj = ufrAll.copy();
ufrAllObj.data = mean(ufrAll.data,2);


figure,
plot(mean(ufrInt(:,[5,6,10]),2));
hold('on');
plot(phz(:,2)*10+50);
plot(uiphz(:)*10+50);


ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'',unitSubset,0.04,true,'gauss');

xyn = xyz.copy();
xyn.data = sq(xyz(:,'nose',[1,2]));
xyl = xyn.copy();;
xyl.resample(lfp);

specArgs = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,40]);
   
[ys,fs,ts] = fet_spec(Trial,lfp,[],[],[],specArgs);

lys = ys.copy;
lys.data = mean(log10(ys(:,fs<15,2)),2);
lys.resample(xyz);

lrs = ys.copy;
lrs.data = mean(log10(ys(:,5<fs&fs<12,3)),2)./mean(log10(ys(:,5<fs&fs<12,4)),2);
lrs.resample(xyz);

uind = sum(ufr.data-eps,2)~=0&nniz(xyz);

tE = decode_bayesian_poisson(smap,(ufr(:,:)'));%,'bins',bins,'baseline',0.001);

%tE(:,:,uind) = tE;
%tE(:,:,~uind) = 0;



% ESTIMATE positions from posterior
epos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
binsE = {};
for i = 1:numel(pfs.adata.bins),  binsE{i} = pfsBins{i}(grows{i});  end
gbins = cell([1,numel(pfsBins)]); 
[gbins{:}] = ndgrid(binsE{:});
gbins = cat(numel(gbins)+1,gbins{:});
ss = substruct('()',[repmat({':'},[1,ndims(gbins)-1]),{1}]);

for tind = 1:size(tE,ndims(tE)),
    ss.subs{end} = tind;
    epos(tind,:) = sum(reshape(gbins.*repmat(subsref(tE,ss),[ones([1,ndims(gbins)-1]),...
                        size(gbins,ndims(gbins))]),...
                               [],size(gbins,ndims(gbins))));
end

% $$$ tE = RectFilter(tE,3,1);
% $$$ tE = permute(RectFilter(permute(tE,[2,1,3]),3,1),[2,1,3]);

tpos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
apos = nan([size(tE,ndims(tE)),1]);
for tind = 1:size(tE,ndims(tE)),
    tbin = LocalMinima2(-tE(:,:,tind),0,30);
    if ~isempty(tbin),
        tpos(tind,:) = [gpfsBins{1}(tbin(1)),gpfsBins{2}(tbin(2))];
        apos(tind)   = tE(tbin(1),tbin(2),tind);
    end
end



stcm = stc2mat(stc,xyz,states);

% $$$ timeSpec = [1:size(ys,1)]./ys.sampleRate;
% $$$ cmapLims = [1,2.5;1,2.5;1,3;1,3.5];

figure();

sp = [];
for s = 1:4,
    sp(s+1) = subplot(6,1,s+1);
    imagesc(ts,fs,log10(ys(:,:,s))');
    axis('xy');
    caxis(cmapLims(s,:));
    colormap(gca,'jet');
end
sp(6) = subplot(6,1,6);
haxSTS = plotSTC(Trial.stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');

linkaxes(get(gcf,'Children'),'x');

sp(1) = subplot(611);
hold('on');
position = animatedline();
position.Marker = '*';
position.MarkerSize = 20;
position.MarkerFaceColor  = 'm';
position.MarkerEdgeColor  = 'm';

predicted = animatedline();
predicted.Marker = '*';
predicted.MarkerSize = 20;
predicted.MarkerFaceColor  = 'm';
predicted.MarkerEdgeColor  = 'm';

for i = 5266163:2:5286163,
    if apos(i)<0.001, continue;end
    
    posterior = imagesc(pfsBins{1}(3:end-2),...
                pfsBins{2}(3:end-2),...
                tE(:,:,i)');
    
    position.clearpoints();
    position.addpoints(xyn(i,1),xyn(i,2),10);
    predicted.clearpoints();
    predicted.addpoints(tpos(i,1),tpos(i,2),10);
    daspect([1,1,1]);
    
    t = NearestNeighbour(timeSpec,stper(i));

    xlim(sp(2),[t-5,t+5]);
    drawnow();
    
    delete(posterior);
end


unitInclusion = sum(logical(ufr(:,:)-eps),2);
statePeriods = any(stcm==7|stcm==7,2);
statePeriods = any(stcm==6|stcm==6,2);

%positionError = sqrt(sum((tpos-xyls).^2,2));
positionError = sqrt(sum((epos-sq(xyz(:,'nose',[1,2]))).^2,2));
lowFreqMeanPower =lys(:);
thetaRatioRadLM =lrs(:);




ind = unitInclusion>=1&round(tpos(:,1)/5)~=0&round(tpos(:,2)/5)~=0&statePeriods&nniz(tpos)&apos>0.001;
figure,
subplot(121);
hist2([lowFreqMeanPower(ind),positionError(ind)],linspace(1,2.5,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])
subplot(122);
hist2([thetaRatioRadLM(ind),positionError(ind)],linspace(0.6,1.1,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])



%% segmentation of immobility sequences ------------------------------------------------------------

mufr = unity(mean(ufr.data,2));



hfig = figure();
hfig.CurrentCharacter = 's';
sp = gobjects([1,2]);;
sp(1) = subplot(221);
plot([1:numel(mufr)]./xyz.sampleRate,mufr);
hold('on');
plot([1:numel(mufr)]./xyz.sampleRate,apos);
plot([1:numel(mufr)]./xyz.sampleRate,unity(ufrAllObj.data));
sp(2) = subplot(223);
haxSTS = plotSTC(stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');


sp(3) = subplot(222);
while ~strcmp(hfig.CurrentCharacter,'q')
    if diff(sp(1).XLim).*xyz.sampleRate<10000,
        cla(sp(3));
        segInd = round(sp(1).XLim.*xyz.sampleRate);
        segInd(segInd<1) = 1;
        segInd(segInd>numel(apos)) = numel(apos);    
        segInd = segInd(1):segInd(2);
        segInd = segInd(apos(segInd)>0.1);
        if ~isempty(segInd),
            axes(sp(3));
            imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
            hold('on');
            plot(tpos(segInd,1),tpos(segInd,2),'-g');
            plot(epos(segInd,1),epos(segInd,2),'-c');            
            scatter(tpos(segInd(1),1),tpos(segInd(1),2),10,'m','filled');
            plot(xyz(segInd(1):segInd(end),'nose',1),xyz(segInd(1):segInd(end),'nose',2),'.','MarkerSize',10);
            axis('xy');
        end
        
        drawnow();
        waitforbuttonpress();
    else
        waitforbuttonpress();
        continue    
    end
end


% Pfeiffer, Foster 2013
% Criteria for spw decoding events
% normalized mean Population must > 1 Hz
% no sequence may jump more than 10cm
%
% output vars
% mean observed speed of each segment
% mean speed of segment


tper = [stc{'theta',xyz.sampleRate}];
%sequencePeriods = ThreshCross(unitInclusion,1,5);
sequencePeriods = [stc{'ripple',xyz.sampleRate}.data];

% REMOVE theta periods
sequencePeriods(WithinRanges(mean(sequencePeriods,2),tper.data),:) = [];


seq.sampleRate = xyz.sampleRate;
seq.periods = sequencePeriods;

seq.duration = nan([size(seq.periods,1),1]);
seq.obsDistance = nan([size(seq.periods,1),1]);
seq.decDistance = nan([size(seq.periods,1),1]);
seq.obsDiffMean = nan([size(seq.periods,1),1]);
seq.decDiffMean = nan([size(seq.periods,1),1]);
seq.obsDiffStd = nan([size(seq.periods,1),1]);
seq.decDiffStd = nan([size(seq.periods,1),1]);
seq.obsDiffMax = nan([size(seq.periods,1),1]);
seq.decDiffMax = nan([size(seq.periods,1),1]);
seq.meanError = nan([size(seq.periods,1),1]);

seq.obsPathDistMean = nan([size(seq.periods,1),1]);
seq.decPathDistMean = nan([size(seq.periods,1),1]);
seq.obsPathDistStd = nan([size(seq.periods,1),1]);
seq.decPathDistStd = nan([size(seq.periods,1),1]);
seq.obsPathAngMean = nan([size(seq.periods,1),1]);
seq.decPathAngMean = nan([size(seq.periods,1),1]);
seq.obsPathAngStd = nan([size(seq.periods,1),1]);
seq.decPathAngStd = nan([size(seq.periods,1),1]);

seq.obsPathAngPPC = nan([size(seq.periods,1),1]);
seq.decPathAngPPC = nan([size(seq.periods,1),1]);


posteriorMaxThr = 0.00001;

dvxy = vxy.data;
dxyz = xyz(:,'nose',[1,2]);
depos = permute(posEstCom,[1,3,2]);
for p = 1:size(seq.periods,1),
    segInd = seq.periods(p,1):seq.periods(p,2);
    segInd = segInd(posteriorMax(segInd)>posteriorMaxThr);

    if numel(segInd)>5,
        seq.duration(p) = diff(segInd([1,end]));
        seq.obsDistance(p) = sqrt(sum(diff(dxyz(segInd([1,end]),1,:)).^2,3));
        seq.decDistance(p) = sqrt(sum(diff(depos(segInd([1,end]),1,:)).^2,3));

        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdistXYZ = sqrt(sum(mdiffXYZ.^2,3));

        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdistEpos = sqrt(sum(mdiffEpos.^2,3));    

        seq.obsPathDistMean(p) = mean(diag(mdistXYZ,1));
        seq.decPathDistMean(p) = mean(diag(mdistEpos,1));

        seq.obsPathDistStd(p) = std(diag(mdistXYZ,1));
        seq.decPathDistStd(p) = std(diag(mdistEpos,1));
        

        obsAng = atan2(diag(mdiffXYZ(:,:,1),1),diag(mdiffXYZ(:,:,2),1));
        decAng = atan2(diag(mdiffEpos(:,:,1),1),diag(mdiffEpos(:,:,2),1));        
        
        seq.obsPathAngMean(p) = circ_mean(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngMean(p) = circ_mean(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        seq.obsPathAngStd(p) = circ_std(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngStd(p) = circ_std(circ_dist(decAng(1:end-1),decAng(2:end)));

        seq.obsPathAngPPC(p) = PPC(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngPPC(p) = PPC(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        
        seq.obsDiffMean(p) = mean(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffMean(p) = mean(mdistEpos(mdistEpos(:)~=0));
        
        seq.obsDiffStd(p) = std(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffStd(p) = std(mdistEpos(mdistEpos(:)~=0));

        if sum(mdistXYZ(:)~=0)~=0,
            seq.obsDiffMax(p) = max(mdistXYZ(mdistXYZ(:)~=0  ));
        else
            seq.obsDiffMax(p) = 0;
        end
        seq.decDiffMax(p) = max(mdistEpos(mdistEpos(:)~=0));
        
        seq.meanError(p)   = mean(sqrt(sum([depos(segInd,1,[1,2])-dxyz(segInd,1,:)].^2,3)));
    end
    
end


figure,
figure();plot(seq.obsPathAngPPC,log10(seq.obsDiffMean),'.')

figure();
hold('on');
plot(seq.obsPathAngPPC,log10(seq.obsPathDistMean),'.')
plot(seq.decPathAngPPC,log10(seq.decPathDistMean),'.')

figure();
subplot(121);
sp = plot(seq.obsPathAngPPC,seq.obsPathAngMean,'.');



subplot(122);
plot(seq.decPathAngPPC,seq.decPathAngMean,'.')



hfig = figure();
sp = gobjects([0,1]);
sp(end+1) = subplot(221);
hold('on');
scatter(seq.decPathAngPPC,log10(seq.decPathDistMean),8,seq.duration,'filled')
sp(end+1) = subplot(222);
hold('on');
scatter(seq.obsPathAngPPC,log10(seq.obsPathDistMean),8,seq.duration,'filled')
linkaxes(get(gcf,'Children'),'xy');
drawnow();
c = gobjects([1,numel(sp)]);;
xy = [0,0];
while isempty(hfig.CurrentCharacter)||hfig.CurrentCharacter~='q',
    waitforbuttonpress();
% REMOVE old unit marker
    delete(c);     
% GET current subplot index    
    axind = find(arrayfun(@isequal,sp,repmat([gca],size(sp))));
% GET xy position of currsor on left mouse button down within current subplot    
    xy = sp(axind).CurrentPoint(1,1:2);
% GET axes Data
    axData = [sp(axind).Children.XData',sp(axind).Children.YData'];
% FIND closest point to currsor on left mouse button down
    [~,mind]=min(sqrt(sum(bsxfun(@minus,xy,axData).^2,2)));
% PLOT Posterior
    segInd = seq.periods(mind,1):seq.periods(mind,2);
    segInd = segInd(posteriorMax(segInd)>0.00001);
    
    subplot(223);
    cla();
    hold('on');
    plot(posEstCom(segInd,1),posEstCom(segInd,2))
    plot(posEstMax(segInd,1),posEstMax(segInd,2))    
    %imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
    %caxis([0,0.75]);
    %colormap('hot');
    axis('xy');
    hold('on');    
    plot(posEstCom(segInd,1),posEstCom(segInd,2),'m','LineWidth',2);
    plot(posEstCom(segInd(1),1),posEstCom(segInd(1),2),'c*','MarkerSize',10);    
    xlim([-500,500]);    
    ylim([-500,500]);        
% $$$     plot(xyn(segInd,1), xyn(segInd,2),'g*','MarkerSize',3);
% $$$     quiver(xyn(segInd,1), xyn(segInd,2),...
% $$$              cos(ang(segInd,'head_back','head_front',1))*100,...
% $$$              sin(ang(segInd,'head_back','head_front',1))*100,0,'Color',[1,1,1]);

    subplot(224);
    cla();    
    hold('on');
    plot(posEstCom(segInd,3),posEstCom(segInd,4));
    plot(posEstMax(segInd,3),posEstMax(segInd,4));
    xlim([-1,1.5]);
    ylim([-0.8,1.8]);
    
    axes(sp(axind));
% HIGHLIGHT selected unit on current subplot
    c(axind) = circle(axData(mind,1),axData(mind,2),0.25,'g');
    for a = find(~ismember(1:numel(sp),axind))
        axes(sp(a));        
        c(a) = circle(sp(a).Children.XData(mind),sp(a).Children.YData(mind),0.25,'k');
    end
end




figure();
plot(seq.obsPathAngPPC,log10(seq.decDiff),'.')
plot(seq.decPathAngPPC,log10(seq.decDiffMax),'.')

figure,
hist2(log10(abs([seq.decMeanDiff,seq.obsDistance])+1e-4),100,100)
Lines([],log10(200),'m');

% Sequence analysis 

% Types of sequences
% THETA 
% PROSPECTIVE 
% SPW 
% 
% Theta Sequence: 
% The hippocampus represents the environment as a squence of action potentials of place cells arranged
% temoporally by the spatial order of the placefields with respect to the rats trajectory. Each sequence
% can last between 
% 
% Filter the sequences based on continuity 





% PROJECTION 

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
ny = cross(nz,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nx = cross(ny,nz);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));

if theta ~= 0,
    evec = cat(2,nx,ny,nz);
    j =1:3;
    headNorm = bsxfun(@rdivide,sq(evec(:,rotationAxisInd,:)),sqrt(sum(evec(:,rotationAxisInd,:).^2,3)));
    headKron = reshape(repmat(headNorm',3,1).*headNorm(:,j(ones(3,1),:)).',[3,3,size(headNorm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    headCPM = reshape(headNorm(:,k)',3,3,size(headNorm,1)).*repmat(j,[1,1,size(headNorm,1)]);

% CREATE rotation matrix
    j = 1:3;
    headRotMat = cos(theta)*repmat(eye(3),[1,1,size(headNorm,1)])...
        +sin(theta)*headCPM...
        +(1-cos(theta))*headKron;

% SET matrix
    ovec = evec(:,1,:);
% ROTATE Basis
    nx = permute(sum(headRotMat.*permute(reshape(nx(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
    ny = permute(sum(headRotMat.*permute(reshape(ny(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
    nz = permute(sum(headRotMat.*permute(reshape(nz(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
end



% DIAGNOSTIC figures ------------------------------------------------------------------------------

% $$$ figure,
% $$$ for u = 1:numel(pfs.data.clu),
% $$$ clf();
% $$$ srmap = pfs.plot(unitSubset(u),'mean',false,[],false,0.25,false,interpParPfsNdim,@jet,mazeMask);
% $$$ subplot2(8,6,1,[1,2]);
% $$$ plot(pft,unitSubset(u),'mean',true,[],false);
% $$$ title(num2str(u));
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ for i = 0:5,
% $$$ for j = 0:6,    
% $$$     subplot2(8,6,j+2,i+1);
% $$$     imagescnan({pfsBins{1},pfsBins{2},srmap(:,:,i*3+20,j*3+15)'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ end
% $$$ end
% $$$ waitforbuttonpress();
% $$$ end
% $$$ mrt = pft.maxRate(unitSubset);



% TEST N-dimensional interpolation and masking of rate map
% $$$ 
% $$$ %srmap = srmap.*mask;
% $$$ figure,
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ subplot(121);
% $$$ imagescnan({pfsBins{1},pfsBins{2},sq(srmap(:,:,25,25))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ subplot(122);
% $$$ imagescnan({pfsBins{3},pfsBins{4},sq(srmap(30,20,:,:))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');

% bhv_decode testing figures

figure,
plot((1:size(fet,1))./fet.sampleRate,fet(:,3));
hold('on');
ind = posteriorMax>0.001;
pt = (1:size(posEstSax,1));
plot(pt(ind)./sampleRate,posEstSax(ind,3),'.');

