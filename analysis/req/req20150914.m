function req20150914(Trial,varargin)
%% tSNE stuff
%mkdir('/storage/gravio/figures/req/req20150914');

initial_dims = 5; % number of dimensions of the feature matrix used
                  % in preliminary rounds of t-SNE

perplexity = 100; % arbitrary scale use to control the malability
                  % of the distributions during dimension reduction

SAVEFIG = true;   % boolean flag, Save figures along with png's 

fetSampleRate = 12;         % Reduce the sample rate of the feature matrix

featureSet = 'fet_tsne_rev5'; % name of the function which constructs
                         % the feature set

ifNormalize = true;     % Flag to normalize using the
                         % z-tranformation

hostPath = '/storage/gravio/figures/'; % Where reportfig should
                                       % save stuff


Trial = MTATrial('jg05-20120317');
Stc = Trial.stc.load(Trial,'hand_labeled_rev3');

[fet,fett,fetd] = feval(featureSet,Trial,fetSampleRate,ifNormalize); % Load Feature matrix of the session


[asmat,labels,keys] =  stc2mat(Stc,fet);                     % Create NxK state matrix
asmat = MTADxyz('data',asmat,'sampleRate',fet.sampleRate);   % Wrap state matrix in MTADxyz object
[~,asmat.data] = max(asmat.data,[],2);                       % Eliminate overlaps by retriving winning state index for each time point
c = jet(numel(Stc.states));                                  % Create base colors
csmat = asmat.copy;                                          % Copy object 
csmat.data = c(csmat.data,:);                                % Fill N time points with corect colors

ind = Trial.stc{'a'}.cast('TimeSeries',fet);                 % Set selection index to all good periods {'gper' a.k.a. 'a'}

start = 1;                                                   % Define start of secondary subselection
stop = size(fet,1);                                          % Define end of secondary subselection 
tind = false([size(fet,1),1]);
tind(start:stop) = true;

skip = 2;                                                    % Third uniform subselection
sind = cat(1,...
           reshape([true(1,round(fet.size(1)/skip));...
                    false(skip-1,round(fet.size(1)/skip))],...
                  [],1),...
           false([mod(fet.size(1),skip),1]));

ind.data = ind.data&sind&tind;                               % Subselection
ind = start:skip:stop;

no_dims = 2;                                                 % N dims for the dimensionallity reduction 
mappedX = tsne(fet(ind,:), csmat(ind,:), no_dims, initial_dims, perplexity);  % Do tSNE

figparm = ['tSNE-fsr_' num2str(fetSampleRate) '-ind_' num2str(start) '_' ...   % Rediculously long name
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-initial_dims_' num2str(initial_dims) ...
            '-no_dims_' num2str(no_dims)];


osts = numel(Stc.states);
hfig = figure(3923923);clf
hold on;
sts = Stc.list_state_attrib('label');
mc = csmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end
legend(sts(1:osts-1));
xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);

reportfig([], hfig, 'tsne', 'req', false,Trial.filebase,figparm,[],SAVEFIG,'png');




%% knn estimated fet maps

mtfet =  MTADxyz('data',mfet(ind,:),'sampleRate',fet.sampleRate);
mtpos =  MTADxyz('data',mappedX,'sampleRate',fet.sampleRate);

[meanMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,125;-125,110],@nanmean);
mmap = reshape(meanMap,numel(Bins{1}),numel(Bins{2}),[]);
[stdMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,125;-125,110],@nanstd);
smap = reshape(stdMap,numel(Bins{1}),numel(Bins{2}),[]);


hfig = figure(38381);
for i = 1:size(mmap,3);
    clf
    [ha,hc] = imagescnan({Bins{1},Bins{2},mmap(:,:,i)},[],false, ...
                        true,[0,0,0]);
    ylabel(hc,'mean z-score','FontName','Courier','FontSize',12);
    axis xy
    title(fett{i},'FontName','Courier','FontSize',12)
    daspect([1,1,1])
   
    reportfig([], hfig, 'tsne', 'req',false,fett{i},fetd{i},[],SAVEFIG,'png',8,8);
end

% SNR maps
% $$$ hfig = figure(38381);
% $$$ for i = 1:size(mmap,3);
% $$$     clf
% $$$     [ha,hc] = imagescnan({Bins{1},Bins{2},mmap(:,:,i)./smap(:,:,i)},[-10,10],false, ...
% $$$                         true,[0,0,0]);
% $$$     ylabel(hc,'mean z-score','FontName','Courier','FontSize',12);
% $$$     axis xy
% $$$     title(fett{i},'FontName','Courier','FontSize',12)
% $$$     daspect([1,1,1])
% $$$    
% $$$     %reportfig([], hfig, 'tsne', 'req',false,fett{i},fetd{i},[],SAVEFIG,'png',8,8);
% $$$ end



