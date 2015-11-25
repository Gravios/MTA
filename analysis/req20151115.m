%function req20151009(Trial)




%% A1 #DESC
% JPDF of the XY distance between the lower spine marker and the
% upper spine marker with contours representing the occupancy of
% other sessions whos labels are derived from the 
%
sltag = 'req20151115';
figTitle = 'JPDF';
figPath  = '/storage/gravio/figures/';
%figPath  = '/gpfs01/sirota/home/gravio/figures/';
baseFigNum = 2838293;

Trials = SessionList(sltag);
numOfTrials = numel(Trials);

Scolors = jet(numOfTrials);
nbin = 50;
tlist = {};
state = {'walk','rear','sit','turn','shake','groom'};
fet = {};
stc = {};

for s = 1:numOfTrials
    tlist{s} = Trials{s}{1};                                   % Load trial list
    Trial = MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2});  % Load MTATrial
    if s ==1,
        [fet{s},ftit,fdesc,fmean,fstd] = fet_tsne(Trial,...       % Load features and save means and stds
                                               Trial.xyz.sampleRate,true);
        stc{s} = Trial.load('stc','hand_labeled_rev2');           % Load States, labeled by hand 

    else
        [fet{s}] = fet_tsne(Trial,...                             % Load features with saved means and stds
                         Trial.xyz.sampleRate,true,fmean,fstd);         
        stc{s} = Trial.load('stc','LGR-hand_labeled_rev2-wrnpms');% Load States, labeled with logistic regression 
    end
end

for i = 1:numel(state),
    k = 1;
    for f = 1:fet{1}.size(2),
        for j = f+1:fet{1}.size(2),
            if f==1&&j==2,
                addBreak = true;
            else
                addBreak = false;
            end
            
            hfig = figure(baseFigNum+k);
            for s = 1:numOfTrials,

                if s == 1,
                    % Clear the figure
                    clf
                    % Index for not zero, inf or nan           
                    ind = nniz(fet{s});                                   
                    % Create nbins over percentile 2 and 98 of fet f
                    edgs    = {linspace([prctile(fet{s}(ind,f),[2,98])+[-2,2],nbin])};
                    % Create nbins over percentile 2 and 98 of fet j
                    edgs(2) = {linspace([prctile(fet{s}(ind,j),[2,98])+[-2,2],nbin])};
                    edc = edgs;
                    [edc{:}] = get_histBinCenters(edc);
                    [X,Y] = meshgrid(edc{:});

                    b = fet{s}(ind,[f,j]);
                    out = hist2(b,edgs{1},edgs{2});
                    imagesc(edgs{1},edgs{2},out');                     % Plot JPDF of hand labeled data
                    axis xy
                end
                
                hold on,
                o = hist2(fet{s}(stc{s}{state{i}},[f,j]),edgs{1},edgs{2});
                F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
                o = conv2(o./sum(o(:)),F,'same');
                contour(X,Y,o',[.0004,.0004],'linewidth',2,'Color',Scolors(s,:))

                if s == numOfTrials,               
                    xlabel(ftit{f});
                    ylabel(ftit{j});
                    title(state{i});
                    caxis([0,1000]);
                    legend(tlist,'location','SouthEastOutSide')
                    hfig.Position  = [100   100   900   650];
                    reportfig(figPath,...             Base path where figures are saved
                    hfig,...                TFigure handel 
                    'CJPDF',...             contours overlayed on joint probability distribution
                    'req',...               Save in the req folder
                    false,...
                        ['fet: ' num2str(f) ' fet: ' num2str(j)],... Tag posted below the img thumbnail
                    [figTitle '-f' num2str(f) 'f' num2str(j)],[],false,'png',[],[],[],addBreak);
                end
            end

            close;
            k=k+1;        
        end
    end

end
