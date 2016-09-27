%% Compute pfk permutations during Ed10-20140820

MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;

set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

sessionList = 'Ed10VR';
  trialList = 'Ed10VR_teleport';

OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';


%% Load and preprocesses data

S = get_session_list(sessionList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');


T = get_session_list(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');
T(3).offsets = [15,-90];
T(1) = [];

if overwriteTrials, QuickTrialSetup(T,'overwrite',true); end

Trial = MTATrial.validate(S(1));

% Create or Load trial epochs
trialLabels = {'Coreg1PEVR1','Shift1PEVR2','Coreg2PEVR3','Shift2PEVR4'};
trialKeys =   {'1','2','3','4'};
if overwriteStc, 
    Trial.stc.updateMode('trialStates');
    Trial.stc.states = {};
    Trial = labelBhvBasic(Trial); 
    Trial = labelTheta(Trial,[],32);

    offsets =  [15,-15;...
                15,-90;...
                15,-15;...
                15,-15];

    Stc = Trial.stc.copy;
    for t = 1:4,
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

units =[1,5,7,9,16,18,22,28,29,99,101,104,107,110,122,134,158,168,184,185]';

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


%% Manually select field centers for double placefields
dunits = [22,29,99,104,134];
dims = 1;
cens = 2;
for d = dunits,
    fprintf('double placefield unit: %i\n',d);
% $$$     for s1 = 1:nsts-1
% $$$         for s2 = s1+1:nsts 

% $$$             pcomx = reshape(pfkpshuff{s1,s2,d==units}.patchCOM(1,:,:,2),[],1);
% $$$             pcomy = reshape(pfkpshuff{s1,s2,d==units}.patchCOM(1,:,:,1),[],1);
% $$$ 
    for t =1:3
        pcomx = reshape(pfkpshuff{t,t+1,d==units}.patchCOM(1,:,:,2),[],1);
        pcomy = reshape(pfkpshuff{t,t+1,d==units}.patchCOM(1,:,:,1),[],1);
            pcind = nniz(pcomx)&nniz(pcomy);
            pcomx = pcomx(pcind);    pcomy = pcomy(pcind);

            figure, plot(pcomx,pcomy,'.');
                    xlim([-600,600,]); 
                    ylim([-350,350,]);            

            [cids] = ClusterPP(gcf);delete(gcf);delete(gcf)

            pcomx = pcomx(cids==1);    pcomy = pcomy(cids==1);

            rsind = randperm(numel(pcomx));
            if numel(rsind)>numIter
                rsind = rsind(1:numIter);
            end

            peakPatchCOMperm(t,t+1,:,d==units,2) = zeros;
            peakPatchCOMperm(t,t+1,1:numel(rsind),d==units,2) = pcomx(rsind);
            peakPatchCOMperm(t,t+1,:,d==units,1) = zeros;
            peakPatchCOMperm(t,t+1,1:numel(rsind),d==units,1) = pcomy(rsind);
% $$$             peakPatchCOMperm(s1,s2,:,d==units,2) = zeros;
% $$$             peakPatchCOMperm(s1,s2,1:numel(rsind),d==units,2) = pcomx(rsind);
% $$$             peakPatchCOMperm(s1,s2,:,d==units,1) = zeros;
% $$$             peakPatchCOMperm(s1,s2,1:numel(rsind),d==units,1) = pcomy(rsind);
            
% $$$         end        
    end
end





%% Bootstrapped place fields
Stc = Trial.stc.copy;
nt = numel(T)-1;
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

for t = 1:4
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


dunits = [22,29,99,104,134];
dims = 1;
cens = 2;
for d = dunits,
    fprintf('double placefield unit: %i\n',d);
    for t = 1:nt

        pcomx = reshape(pfkshuff{t,d==units}.patchCOM(1,:,:,2),[],1);
        pcomy = reshape(pfkshuff{t,d==units}.patchCOM(1,:,:,1),[],1);
        pcind = nniz(pcomx)&nniz(pcomy);
        pcomx = pcomx(pcind);
        pcomy = pcomy(pcind);

        figure,
        plot(pcomx,pcomy,'.')
        xlim([-600,600,]);        
        ylim([-350,350,]);            
        [cids] = ClusterPP(gcf);    
        delete(gcf)
        
        pcomx = pcomx(cids==1);
        pcomy = pcomy(cids==1);
        
        rsind = randperm(numel(pcomx));
        if numel(rsind)>numIter
            rsind = rsind(1:numIter);
        end
        
        peakPatchCOM(t,:,d==units,2) = zeros;
        peakPatchCOM(t,1:numel(rsind),d==units,2) = pcomx(rsind);
        peakPatchCOM(t,:,d==units,1) = zeros;
        peakPatchCOM(t,1:numel(rsind),d==units,1) = pcomy(rsind);
        
        delete(gcf)
        
    end
end



% Calculate the dprime for each distribit
dprime = @(x,y) (mean(x(nniz(x')))-mean(y(nniz(y'))))/sqrt(0.5*(var(x(nniz(x')))+var(y(nniz(y')))));
dprx = nan([nt-1,numel(units)]);
dpry = nan([nt-1,numel(units)]);

for i = 1:3,
    for u = 1:numel(units);    
        dprx(i,u) = dprime(peakPatchCOM(i,:,u,2),peakPatchCOM(i+1,:,u,2));
        dpry(i,u) = dprime(peakPatchCOM(i,:,u,1),peakPatchCOM(i+1,:,u,1));
    end    
end




mRate = [];
for t = 1:nt;
    for u = units'
        mRate(t,find(u==units)) = pfk{t}.maxRate(u);
    end
end





FigDir = 'Ed10-20140820-permutation_test_selected2';
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

    for t = 1:4,
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
        
        
        if t<4,
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

    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_perm_selected',num2str(units(u)),'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_perm_selected',num2str(units(u)),'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);
end

save('/storage/gravio/ownCloud/Shared/VR_Methods/matlab/pfs_shift_perm_0820_all.mat','-v7.3');




save('/storage/gravio/ownCloud/Shared/VR_Methods/matlab/pfs_shift_perm_0820.mat','-v7.3',...
     'units',...
     'peakPatchCOM',...
     'peakPatchRate',...
     'permPval',...
     'dprx',...
     'dpry');
