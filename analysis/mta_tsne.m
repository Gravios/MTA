function out = mta_tsne(Trial,varargin)
% function mta_tsne(Trial,varargin)
%
%    varargin:
%
%           Fet ('fet_tsne_rev15')
%           sampleRate (12)    
%           Stc ('hand_labeled')    
%           states({'rear','walk','turn','pause','groom','sit'})
%           initDims (5)    
%           nDims (2)
%           perplexity (100)    
%           ifNorm (true)
%           ifReportFig (false)
%           subset ([])
%           hostPath ('/storage/gravio/figures/')
%           overwrite (false)
%

defArgs = {                                                 ...
    ...        Fet
              'fet_tsne_rev15',                             ...
    ...
    ...        sampleRate
               12,                                          ...
    ...
    ...        Stc
               'hand_labeled',                              ...
    ...
    ...        states
               {'rear','walk','turn','pause','groom','sit'},...
    ...
    ...        initDims
               5,                                           ...
    ...
    ...        nDims
               2,                                           ...
    ...
    ...        perplexity
               100,                                         ...
    ...
    ...        ifNorm
               true,                                        ...
    ...
    ...        ifReportFig
               false,                                       ...
    ...
    ...        subset
               [],                                          ...
    ...
    ...        hostPath
               '/storage/gravio/figures/',                  ...
    ...
    ...        overwrite
               false                                        ...
};


[Fet,sampleRate,Stc,states,initDims,nDims,...
 perplexity,ifNorm,ifReportFig,subset,hostPath,overwrite] = DefaultArgs(varargin,defArgs);


% Load Features
if ischar(Fet),
    [Fet,fett,fetd] = feval(fet,Trial,sampleRate,ifNorm); % Load Feature matrix of the session
elseif isa(Fet,'MTADfet')
    sampleRate = Fet.sampleRate; 
else
    error('MTA:analysis:mta_tsne:UnknownFet');    
end

% Load State Collection
if ischar(Stc)
    Stc = Trial.load('stc',Stc);
elseif ~isa(Stc,'MTAStateCollection'),
    error('MTA:analysis:mta_tsne:UnknownStc');    
end
    

% Create Color Mapping 
[asmat,labels,keys] =  stc2mat(Stc,Fet,states);              % Create NxK state matrix
asmat = MTADxyz('data',asmat,'sampleRate',Fet.sampleRate);   % Wrap state matrix in MTADxyz object
[~,asmat.data] = max(asmat.data,[],2);                       % Eliminate overlaps by retriving winning state index for each time point
c = jet(numel(Stc.states));                                  % Create base colors
csmat = asmat.copy;                                          % Copy object 
csmat.data = c(csmat.data,:);                                % Fill N time points with corect colors


% Select Range
ind = Stc{'a'}.cast('TimeSeries');                           % Set selection index - good periods {'gper','a'}
ind.resample(Fet);

if isempty(subset)                                           % Second subselection
    subset = [1,Fet.size(1)];
end
tind = false([size(Fet,1),1]);
tind(subset(1):subset(2)) = true;

skip = 2;                                                    % Third uniform subselection
sind = cat(1,...
           reshape([true(1,round(Fet.size(1)/skip));...
                    false(skip-1,round(Fet.size(1)/skip))],...
                  [],1),...
           false([mod(Fet.size(1),skip),1]));

ind.data = ind.data&sind&tind;                               % Subselection


% Build filepath
figparm = ['tSNE-' Fet.label...
           '_sr_' num2str(sampleRate) ...
           '_subset_' num2str(subset(1)) '-' num2str(skip) '-' num2str(subset(2)) ...
           '_perplexity_' num2str(perplexity) ...
           '_initDims_' num2str(initDims) ...
           '_nDims_' num2str(nDims)];
filepath = fullfile(Trial.path.data,'analysis',figparm);


% Run tSNE or Load data
if ~exist(filepath,'file')||overwrite
    if ifReportFig, % lazy
        mappedX = tsne(Fet(ind,:), csmat(ind,:), nDims, initDims, perplexity); 
    else
        mappedX = tsne(Fet(ind,:), [], nDims, initDims, perplexity);
    end
    save(fullfile(Trial.path.data,'analysis',figparm),'states','Fet','csmat','ind','mappedX','perplexity');
else
    load(filepath);
end


% Don't complain about this I won't hear it
if nargout>0
    out = load(filepath);
end

% Report Figure 
if ifReportFig, 
    osts = numel(states);
    hfig = figure(3923923);clf
    hold on;
    mc = csmat(ind,:);
    for nc = 1:osts,
        nind = all(bsxfun(@eq,c(nc,:),mc),2);
        h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
        try,h.MarkerFaceColor = h.CData(1,:);end
    end
    legend(states);
    xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
    ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);

    reportfig([], hfig, 'tsne', 'req', false,Trial.filebase,figparm,[],true,'png');
end


