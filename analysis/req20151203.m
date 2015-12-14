function req20151203(Trial,varargin)
%function req20151203(Trial,varargin)
% tSNE over all jg05 sessions excluding jg05-20120317
% jg05-20120317 features are then mapped onto tsne space

%mkdir('/storage/gravio/figures/req/req20151203');


initial_dims = 5;       % number of dimensions of the feature matrix used
                        % in preliminary rounds of t-SNE
perplexity = 100;       % arbitrary scale use to control the malability
                        % of the distributions during dimension reduction
SAVEFIG = true;         % boolean flag, Save figures along with png's 
NEW_SAMPLE_RATE = 10;   % Reduce the sample rate of the feature matrix
HOST_PATH = '/storage/gravio/figures/'; % Where reportfig should
                                       % save stuff


%Trial = MTATrial('jg05-20120317');
%Trial.stc.load(Trial,'hand_labeled_rev2');
%Trial = MTATrial('Ed01-20140707');
%Stc = Trial.stc.load(Trial,'hand_labeled_rev1');
slist = SessionList('req20151203');


%Reference Trial Stuff
RefTrial = MTATrial('jg05-20120317');
RefState = RefTrial.stc{'a'};
rfet = fet_tsne(RefTrial,NEW_SAMPLE_RATE);

for s = 1:numel(slist),
    Trial = MTATrial(slist{s}{1},slist{s}{3},slist{s}{2});

  % Load Feature matrix of the session    
  %%% Fix unity 
    [tfet,fett,fetd] = fet_tsne(Trial,NEW_SAMPLE_RATE);
    [~,Rmean,Rstd] = unity(rfet,[],[],[],[]);
    tfet.normalize_to_reference(rfet);
    tfet.unity(rfet,[],Rmean,Rstd);

    if s == 1,
        fet = tfet.copy;
    else
        fet.data = cat(1,fet.data,tfet(Trial.stc{'a'},:));
    end
    
end
fet.data(~nniz(fet),:) = [];

%Stc = Trial.stc.load(Trial,'hand_labeled_rev2');





% $$$ [asmat,labels,keys] =  stc2mat(Stc,fet);
% $$$ asmat = MTADxyz('data',asmat,'sampleRate',fet.sampleRate);
% $$$ [~,asmat.data] = max(asmat.data,[],2);
% $$$ c = jet(numel(Stc.states));
% $$$ csmat = asmat.copy; 
% $$$ csmat.data = c(csmat.data,:);
% $$$ 
% $$$ ind = Trial.stc{'a'};
% $$$ 
% $$$ mfet = fet(ind,:);
% $$$ msmat = csmat(ind,:);

start = 1;
skip = 5;
stop = size(fet,1);
no_dims = 2;

ind = start:skip:stop;
mappedX = tsne(fet(ind,:), [], no_dims, initial_dims, perplexity); 



%pos = map_feature_to_tsne_space(Fet,tsneFet,tsneMap)

%%%% Save this next time
comptSneMap = zeros([fet.size(1),2]);
tic
for s = 1:fet.size(1),
    [mdist,mind] = sort(sqrt(sum(bsxfun(@minus,fet(ind,:),fet(s,:)).^2,2)));
    dwght = 1./(mdist(1:10)+eps)./sum(1./(mdist(1:10)+eps));
    comptSneMap(s,:) = sum(mappedX(mind(1:10),:).*repmat(dwght,[1,size(mappedX,2)]));
    if ~mod(s,10000),toc,disp([num2str(s) ' of ' num2str(fet.size(1))]),tic,end
end
toc



%%Start Here





figparm = ['tSNE-nsr_' num2str(NEW_SAMPLE_RATE) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-initial_dims_' num2str(initial_dims) ...
            '-no_dims_' num2str(no_dims)];




osts = numel(Stc.states);
hfig = figure(3923923);clf
hold on;
sts = Stc.list_state_attrib('label');
mc = msmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end
legend(sts(1:osts-1));
xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);

reportfig([], hfig, 'tsne', 'req', false,Trial.filebase,figparm,[],SAVEFIG,'png');


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
