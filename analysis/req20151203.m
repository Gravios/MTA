%function req20151203(Trial,varargin)
%function req20151203(Trial,varargin)
% tSNE over all jg05 sessions excluding jg05-20120317
% jg05-20120317 features are then mapped onto tsne space

%mkdir('/storage/gravio/figures/req/req20151203');


initial_dims = 3;       % number of dimensions of the feature matrix used
                        % in preliminary rounds of t-SNE
perplexity = 80;       % arbitrary scale use to control the malability
                        % of the distributions during dimension reduction
SAVEFIG = false;         % boolean flag, Save figures along with png's 
NEW_SAMPLE_RATE = 10;   % Reduce the sample rate of the feature matrix
HOST_PATH = '/storage/gravio/figures/'; % Where reportfig should
                                       % save stuff
states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_tsne_rev15';
normalize = true;
map2reference = true;
sessionSet = 'hand_labeled';
mfilename = 'req20151203';


%Reference Trial Stuff
RefTrial = MTATrial('jg05-20120317');
if normalize,
    RefState = RefTrial.stc{'a'};
    rfet = feval(featureSet,RefTrial,NEW_SAMPLE_RATE);
    rfet.data = rfet(RefState,:);
    [~,Rmean,Rstd] = unity(rfet,[],[],[],[]);
end


if map2reference, mapping =    ['-map2_' RefTrial.filebase];else mapping    = '';end
if normalize,     normStatus = '-norm';                     else normStatus = '';end

fileLoc = fullfile(MTASession([]).path.data,'analysis',...
                   [mfilename,'-',sessionSet,'-',featureSet,mapping,normStatus,'.mat']);

if ~exist(fileLoc,'file'),

    slist = SessionList(sessionSet);

    sfet = [];
    Stc = {};
    sts = [];
    for s = 1:numel(slist),
        Trial = MTATrial(slist(s).sessionName,slist(s).trialName,slist(s).mazeName);
        Trial.load('stc',slist(s).stcMode);
        Stc = Trial.stc.copy;
        % Load Feature matrix of the session    
        [tfet,fett,fetd] = feval(featureSet,Trial,NEW_SAMPLE_RATE);

        if strcmp(Trial.filebase,RefTrial.filebase)&&map2reference,
            tfet.map_to_reference_session(Trial,RefTrial);
        end

        if s == 1,
            fet = tfet.copy;
            fet.data = fet(Stc{'a'},:);
        else
            fet.data = cat(1,fet.data,tfet(Stc{'a'},:));
        end

        
        tsts = [];
        for t = states,
            tper = Stc{t{1}};
            tper.cast('TimeSeries');
            tper.resample(fet.sampleRate);
            tsts = cat(2,tsts,tper(Stc{'a'}));
        end     
        % Inefficient ... I don't care
        sts = cat(1,sts,tsts);
        sts(~nniz(fet),:) = [];        
        fet.data(~nniz(fet),:) = [];
        sfet(end+1) = length(fet.data);
    end



    if normalize
        fet.unity([],Rmean,Rstd);
    end

    

    % Add colors to states
    asmat = MTADfet('data',sts,'sampleRate',fet.sampleRate);
    [~,asmat.data] = max(asmat.data,[],2);
    c = jet(numel(states));
    c = [0,0,1;...
         1,0,0;...
         0,1,0;...
         0,1,1;...
         1,0,1;...
         1,1,0;];
    csmat = asmat.copy; 
    csmat.data = c(csmat.data,:);





    start = 1;
    skip = 3;
    stop = size(fet,1);
    no_dims = 2;
    ind = start:skip:stop;


    mappedX = tsne(fet(ind,:), csmat(ind,:), no_dims, initial_dims, perplexity); 

    
    save(fileLoc);
else
    load(fileLoc);
end


%% For when all points have state labels

osts = numel(states);
hfig = figure(3923924);clf
hold on;
mc = csmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end
legend(states,'location','south','Orientation','horizontal');
reportfig(HOST_PATH,           ...
          'FigHandle',hfig,    ...
          'FileName','tsne',   ...
          'FigDir','req',      ...
          'Preview',false,     ...
          'Tag',[sessionSet,'-',featureSet,'-nomap'],...
          'Comment','tsne',    ...
          'Resolution',100,    ...
          'SaveFig', SAVEFIG,  ...
          'format','png',      ...
          'width',8,           ...
          'height',8);


