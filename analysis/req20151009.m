function req20151009(Trial)




%% A1 #DESC
% JPDF of the XY distance between the lower spine marker and the
% upper spine marker with contours representing the occupancy of
% other sessions whos labels are derived from the 
%
sltag = 'req20151009';                                         
figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
figPath  = '/storage/gravio/figures/';
figPath  = '/gpfs01/sirota/home/gravio/figures/';
baseFigNum = 2838293;

Trials = SessionList(sltag);
numOfTrials = numel(Trials);

Scolors = jet(numOfTrials);
nbin = 75;
tlist = {};
state = 'rear';

for s = 1:numOfTrials
    tlist{s} = Trials{s}{1};                                   % Load trial list
    Trial = MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2});  % Load MTATrial
    if s ==1,
        [fet,ftit,fdesc,fmean,fstd] = fet_tsne(Trial,...       % Load features and save means and stds
                                               Trial.xyz.sampleRate,true);
        stc = Trial.load('stc','hand_labeled_rev2');           % Load States, labeled by hand 

    else
        [fet] = fet_tsne(Trial,...                             % Load features with saved means and stds
                         Trial.xyz.sampleRate,true,fmean,fstd);         
        stc = Trial.load('stc','LGR-hand_labeled_rev2-wrnpms');% Load States, labeled with logistic regression 
    end
    %req20151009_A1(Trial,hfig,state,Scolors(s,:));             % Plot contours for each feature


    k = 1;
    for f = 1:fet.size(2),
        for j = f+1:fet.size(2),
            hfig = figure(baseFigNum+k);
            k=k+1;

            if s == 1,
                % Clear the figure
                clf
                % Index for not zero, inf or nan           
                ind = nniz(fet);                                   
                % Create nbins over percentile 2 and 98 of fet f
                edgs    = {linspace([prctile(fet(ind,f),[2,98]),nbin])};
                % Create nbins over percentile 2 and 98 of fet j
                edgs(2) = {linspace([prctile(fet(ind,j),[2,98]),nbin])};
                edc = edgs;
                [edc{:}] = get_histBinCenters(edc);
                [X,Y] = meshgrid(edc{:});

                b = fet(ind,[f,j]);
                out = hist2(b,edgs{1},edgs{2});
                imagesc(edgs{1},edgs{2},out');                     % Plot JPDF of hand labeled data
                axis xy
            end

            hold on,
            o = hist2(fet(stc{state},[f,j]),edgs{1},edgs{2});
            F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
            o = conv2(o./sum(o(:)),F,'same');
            contour(X,Y,o',[.0004,.0004],'linewidth',2,'Color',Scolors(s,:))
            
        end
    end

end


for k = 1:fet.size(2),
    hfig = figure(baseFigNum+k);k=k+1;
    caxis([0,200])
    legend(tlist,'location','SouthEast')
    hfig.Position  = [100   100   782   629];
    reportfig(figPath,...             Base path where figures are saved
    hfig,...                Figure handel 
    'CJPDF',...             contours overlayed on joint probability distribution
    'req',...               Save in the req folder
    false,...
        [sltag '_A1'],...
        figTitle,[],false,'png',[],[],[],false);
end

% $$$ 
% $$$ 
% $$$ %% A2 #DESC
% $$$ %
% $$$ %
% $$$ sltag = 'req20151009';
% $$$ figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
% $$$ figPath  = '/storage/gravio/figures/';
% $$$ hfig = figure(2838293);clf
% $$$ Trials = SessionList(sltag);
% $$$ numOfTrials = numel(Trials);
% $$$ Scolors = jet(numOfTrials);
% $$$ tlist = {};
% $$$ state = 'rear';
% $$$ for s = 1:numOfTrials
% $$$     if s == 1, AddBreak = true; else  AddBreak = false; end
% $$$     clf(hfig);
% $$$     tlist{s} = Trials{s}{1};
% $$$     req20151009_A2(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state);
% $$$     axis xy
% $$$     caxis([0,200])
% $$$     legend({state},'location','SouthEast')
% $$$     hfig.Position  = [100   100   782   629];
% $$$     reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A2'],figTitle,[],false,'png',[],[],[],AddBreak);
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ %% B1 #DESC
% $$$ % JPDF of the XY distance between the lower spine marker and the
% $$$ % upper spine marker with contours representing the occupancy of
% $$$ % other sessions whos labels are derived from the 
% $$$ %
% $$$ sltag = 'req20151009';
% $$$ figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
% $$$ figPath  = '/storage/gravio/figures/';
% $$$ hfig = figure(2838293);clf
% $$$ Trials = SessionList(sltag);
% $$$ numOfTrials = numel(Trials);
% $$$ Scolors = jet(numOfTrials);
% $$$ tlist = {};
% $$$ state = 'rear';
% $$$ for s = 1:numOfTrials
% $$$     tlist{s} = Trials{s}{1};
% $$$     req20151009_A1(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state,Scolors(s,:));
% $$$ end
% $$$ caxis([0,200])
% $$$ legend(tlist,'location','SouthEast')
% $$$ hfig.Position  = [100   100   782   629];
% $$$ reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A1'],figTitle,[],false,'png',[],[],[],false);
% $$$ 
% $$$ 
% $$$ %% B2 #DESC
% $$$ %
% $$$ %
% $$$ sltag = 'req20151009';
% $$$ figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
% $$$ figPath  = '/storage/gravio/figures/';
% $$$ hfig = figure(2838293);clf
% $$$ Trials = SessionList(sltag);
% $$$ numOfTrials = numel(Trials);
% $$$ Scolors = jet(numOfTrials);
% $$$ tlist = {};
% $$$ state = 'rear';
% $$$ for s = 1:numOfTrials
% $$$     if s == 1, AddBreak = true; else  AddBreak = false; end
% $$$     clf(hfig);
% $$$     tlist{s} = Trials{s}{1};
% $$$     req20151009_A2(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state);
% $$$     caxis([0,200])
% $$$     legend({state},'location','SouthEast')
% $$$     hfig.Position  = [100   100   782   629];
% $$$     reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A2'],figTitle,[],false,'png',[],[],[],AddBreak);
% $$$ end
% $$$ 
% $$$ 
% $$$ 
