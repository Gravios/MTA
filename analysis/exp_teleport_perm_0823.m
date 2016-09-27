%% Compute pfk permutations during Ed10-20140820

MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;

set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

sessionList = 'Ed10VR';
  trialList = 'Ed10VR_20160823';

OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';


%% Load and preprocesses data

S = SessionList(sessionList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');

T = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');
T(1) = [];

if overwriteTrials, QuickTrialSetup(T,'overwrite',true); end

Trial = MTATrial.validate(S(4));

% Create or Load trial epochs
trialLabels = {'Coreg1PEVR1','Shift1PEVR2','starfield3'};
trialKeys =   {'1','2','3'};
if overwriteStc, 
    Trial.stc.updateMode('trialStates');
    Trial.stc.states = {};
    Trial = labelBhvBasic(Trial); 
    Trial = labelTheta(Trial,[],32);

    offsets =  [15,-15;...
                15,-15;...
                15,-15];

    Stc = Trial.stc.copy;
    for t = 1:numel(trialLabels),
        aper = resample(Trial.sync.copy,Trial.xyz.sampleRate);
        aper.data = aper.data(t,:)-aper.data(1)+1;
        aper = aper+offsets(t,:);
        Stc.addState(Trial.spath,...
                     Trial.filebase,...
                     aper.data,...
                     Trial.xyz.sampleRate,...
                     Trial.sync.copy,...
                     Trial.sync.data(1),...
                     trialLabels{t},...
                     trialKeys{t});
    end
    Stc.updateMode('trialStates');
    Stc.save(1);
    
else
    Stc = Trial.load('stc','trialStates');
end

Trial.stc = Stc;

states = cellfun(@horzcat,...
                 trialLabels,...
                 repmat({'&velHthresh'},size(trialLabels)),...
                 'UniformOutput',false);

units = [13,17,18,20,21,28,32,42,48,56,57,58,59,62,64,68,77,80,87,107,112,116,120,142,147,160,164,167,172,174];


nsts = numel(states);

binDims = [20,20];
numIter = 1001;
nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 2;
sampleRate = 30;
pfk = {};
overwrite = false;

xyz = Trial.load('xyz');
xyz.resample(sampleRate);

for s1 = 1:nsts-1
    for s2 = s1+1:nsts
        pfkp{s1,s2} = MTAAknnpfs_perm(Trial,units,states([s1,s2]),overwrite, ...
                                     'binDims',binDims,...
                                     'nNearestNeighbors',nNearestNeighbors,...
                                     'ufrShufBlockSize',ufrShufBlockSize,...
                                     'distThreshold',distThreshold,...
                                     'pos',xyz,...
                                     'numIter',numIter);
    end
end

        


pfkpstats = {};
pfkpshuff = {};
% Test this version should be able to run multiple units at once
for s1 = 1:nsts-1
    for s2 = s1+1:nsts
        for u = 1:numel(units)        
            [pfkpstats{s1,s2,u},pfkpshuff{s1,s2,u}] = PlaceFieldStats(Trial,pfkp{s1,s2},units(u));
        end
    end
end



peakPatchCOMperm  = nan([nsts,nsts,numIter,numel(units),2]);
peakPatchRatePerm = nan([nsts,nsts,numIter,numel(units)]);
peakPatchAreaPerm = nan([nsts,nsts,numIter,numel(units)]);
for s1 = 1:nsts-1
    for s2 = s1+1:nsts
        for k = 1:numIter,
            % Retrieve the patch center of mass from patch with the highest firing rate
            %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
            pcom = ...
            sq(cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                    pfkpshuff(s1,s2,:),...
                    repmat({k},[1,1,numel(units)]),...
                    'UniformOutput',false))';

            pind = ~cellfun(@isempty,pcom);
            peakPatchCOMperm(s1,s2,k,pind,:) = ...
                cell2mat(cellfun(@(x) x(:,1),...
                                 pcom(pind),...
                                 'uniformoutput',false))';


            
            
            % Retrieve the max patch Firing rate
            peakPatchRatePerm(s1,s2,k,:) = ...
                sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                           pfkpshuff(s1,s2,:),...
                           repmat({k},[1,1,numel(units)])));

            
            % Retrieve the patch area from patch with highest firing rate
            parea = ...
                sq(cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                        pfkpshuff(s1,s2,:),...
                        repmat({k},[1,1,numel(units)]),...
                        'UniformOutput',false))';
            pind = ~cellfun(@isempty,parea);
            peakPatchAreaPerm(s1,s2,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                              parea(pind),...
                                                              'uniformoutput',false))');
        end
    end
end




%% Bootstrapped place fields
Stc = Trial.stc.copy;
nt = numel(T);
states = {'velHthresh'};
nsts = numel(states);

binDims = [20,20];
numIter = 1001;
nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 2;
sampleRate = 30;
pfk = {};
overwrite = false;

for t = 1:nt
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial); 
    xyz = Trial.load('xyz');
    xyz.resample(sampleRate);

    pfk{t} = MTAAknnpfs_bs(Trial,units,states,overwrite, ...
                             'binDims',binDims,...
                             'nNearestNeighbors',nNearestNeighbors,...
                             'ufrShufBlockSize',ufrShufBlockSize,...
                             'distThreshold',distThreshold,...
                             'pos',xyz,...
                             'numIter',numIter);
end

pfkstats = {};
pfkshuff = {};
% Test this version should be able to run multiple units at once
for t = 1:nt
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial);
    for u = 1:numel(units)        
        [pfkstats{t,u},pfkshuff{t,u}] = PlaceFieldStats(Trial,pfk{t},units(u));
    end
end

for t = 1:nt
    for k = 1:numIter,
        % Retrieve the patch center of mass from patch with the
        % highest firing rate
        %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...        
        pcom = ...
            cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                    pfkshuff(t,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);

        pind = ~cellfun(@isempty,pcom);
        peakPatchCOM(t,k,pind,:) = ...
            cell2mat(cellfun(@(x) x(:,1),...
                             pcom(pind),...
                             'uniformoutput',false))';

        
        % Retrieve the max patch Firing rate
        peakPatchRate(t,k,:) = ...
            sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                       pfkshuff(t,:),...
                       repmat({k},[1,numel(units)])));

        
        % Retrieve the patch area from patch with highest firing rate
        parea = ...
            cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                    pfkshuff(t,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);
        pind = ~cellfun(@isempty,parea);
        peakPatchArea(t,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                      parea(pind),...
                                                      'uniformoutput',false))');
    end
end






mRate = [];
for t = 1:nt;
    for u = units
        mRate(t,find(u==units)) = pfk{t}.maxRate(u);
    end
end


% Calculate the dprime for each distribit
dprime = @(x,y) (mean(x(nniz(x')))-mean(y(nniz(y'))))/sqrt(0.5*(var(x(nniz(x')))+var(y(nniz(y')))));
dprx = nan([nt-1,numel(units)]);
dpry = nan([nt-1,numel(units)]);

for i = 1:2,
    for u = 1:numel(units);    
        dprx(i,u) = dprime(peakPatchCOM(i,:,u,2),peakPatchCOM(i+1,:,u,2));
        dpry(i,u) = dprime(peakPatchCOM(i,:,u,1),peakPatchCOM(i+1,:,u,1));
    end    
end



FigDir = 'Ed10-20140823-permutation_test_selected2';
%FigDir = [Trial.filebase,'-permutation_test'];
mkdir(fullfile(OwnDir,FigDir))


figHnum = 399329240;;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[5 0 27 25])
set(hfig,'PaperPositionMode','auto');

autoincr = true;
unit =units(1);

eds = linspace(-1000,1000,500);
while unit~=-1
    clf

    for t = 1:nt,
        ypos = [2*t-1,2*t-1+1];

        % Plot Place Fields
        subplot2(8,2,ypos,1);
        hold('on');
        pf = pfk{t};        
        % Correct color of nans and plot place field
        ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,1),fliplr(pf.adata.binSizes'));
        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
        % Add text: max rate
        text(pf.adata.bins{1}(end)-200,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        caxis([-1,sq(max(mRate(:,unit==units)))]);
        caxis([-1,10]);        
        ylabel([pf.session.trialName ':' pf.parameters.states]);

        plot(median(peakPatchCOM(t,:,unit==units,2)),...
             median(peakPatchCOM(t,:,unit==units,1)),'*k');
        xlim([-600,600]),ylim([-350,350])                    
        if t==1,
            title(['Unit: ',num2str(unit)]);
        end
        
        
        if t<nt,
            % Calculate distributions
            u = find(unit==units);
            s = [t,t+1];

            ind2 = peakPatchCOM(s(2),:,u,2)~=0;
            ind1 = peakPatchCOM(s(1),:,u,2)~=0;

            cpdiff = bsxfun(@minus,peakPatchCOM(s(2),ind2(2:end),u,2),peakPatchCOM(s(1),ind1(2:end),u,2)');
            cpdiff = cpdiff(:);

            gpdiff = bsxfun(@minus,sq(peakPatchCOMperm(s(1),s(2),2:2:end,u,2)),sq(peakPatchCOMperm(s(1),s(2),3:2:end,u,2))');
            gpdiff = gpdiff(:);
            
            %plot permutation distribution 

            sp = subplot2(8,2,ypos+1,2);
            sp.Units = 'centimeters';

            sp.Position(4) = sp.Position(4)-1;
            
            hold('on');
            hs = bar(eds,histc(gpdiff,eds),'histc');
            hs.FaceColor = 'c';
            hs.EdgeColor = 'c';
            hs.FaceAlpha = 0.5;
            hs.EdgeAlpha = 0.5;

            hs = bar(eds,histc(cpdiff,eds),'histc');
            hs.FaceColor = 'm';
            hs.EdgeColor = 'm';
            hs.FaceAlpha = 0.2;
            hs.EdgeAlpha = 0.2;

            Lines(median(cpdiff),[],'g');

            if mod(t,2),
                compFunc = @gt;
            else
                compFunc = @lt;
            end
            
            permPval(t,u) = sum(compFunc(median(cpdiff),gpdiff))/numel(gpdiff);

% $$$             title({['Permutation test ',pfk{t}.session.trialName,...
% $$$                     ' VS ',...
% $$$                     pfk{t+1}.session.trialName],... 
% $$$                    ['pval: ',num2str(permPval(u))]});
            title(['pval: ',num2str(permPval(t,u))]);
            
        end    
    end

    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_perm',num2str(units(u)),'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_perm',num2str(units(u)),'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);
end





figure,plot(dprx(2,:),log(permPval(2,:)),'.'),Lines([],-3,'k');

save('/storage/gravio/ownCloud/Shared/VR_Methods/matlab/pfs_shift_perm_0823_all.mat','-v7.3');


save('/storage/gravio/ownCloud/Shared/VR_Methods/matlab/pfs_shift_perm_0823.mat','-v7.3',...
     'units',...
     'peakPatchCOM',...
     'peakPatchRate',...
     'permPval',...
     'dprx',...
     'dpry');
