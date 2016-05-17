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

defArgs = struct('Fet',           'fet_tsne_rev15',...
                 'sampleRate',    12,...
                 'Stc',           'hand_labeled',...
                 'states',        {{'rear','walk','turn','pause','groom','sit'}},...
                 'initDims',      5,...
                 'nDims',         2,...
                 'perplexity',     100,...
                 'ifNorm',        true,...
                 'ifReportFig',   false,...
                 'subset',        [],...
                 'overwrite',      false);

[Fet,sampleRate,Stc,states,initDims,nDims,...
 perplexity,ifNorm,ifReportFig,subset,overwrite] = DefaultArgs(varargin,defArgs,'--struct');

% Load Features
if ischar(Fet),
    [Fet,fett,fetd] = feval(fet,Trial,sampleRate,ifNorm); % Load Feature matrix of the session
elseif isa(Fet,'MTADfet')
    sampleRate = Fet.sampleRate; 
    if Fet.isempty,
        Fet.load(Trial);
    end
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
aind = sum(asmat,2)~=0;
asmat = MTADfet('data',asmat,'sampleRate',Fet.sampleRate);   % Wrap state matrix in MTADxyz object
[~,asmat.data] = max(asmat.data,[],2);                       % Eliminate overlaps by retriving winning state index for each time point
switch numel(states)
  case 1
    c = [0,0,1];
  case 2
    c = [0,0,1;...
         1,0,0];
  case 3
    c = [0,0,1;...
         0,1,0;...         
         1,0,0];
  otherwise
    c = jet(numel(Stc.states));                                  % Create base colors
end

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

ind.data(isnan(ind.data))=0;
ind.data = ind.data&sind&tind&aind&nniz(Fet);                          % Subselection


% Build filepath
figparm = ['tSNE-' Fet.label...
           '_sr_' num2str(sampleRate) ...
           '_subset_' num2str(subset(1)) '-' num2str(skip) '-' num2str(subset(2)) ...
           '_perplexity_' num2str(perplexity) ...
           '_initDims_' num2str(initDims) ...
           '_nDims_' num2str(nDims),...
           '_states_' strjoin(states,'_')];
filepath = fullfile(Trial.path.data,'analysis',[figparm,'.mat']);


% Run tSNE or Load data
if ~exist(filepath,'file')||overwrite
    mappedX = tsne(Fet(ind,:), [], nDims, initDims, perplexity);
    save(filepath,'states','mappedX','perplexity',...
                  'filepath','ind','csmat','nDims','initDims')
end
load(filepath);

% Don't complain about this I won't hear it
if nargout>0
    out = load(filepath);
end

% Report Figure 
if ifReportFig, 
    folderName = dbstack;
    if numel(folderName)>1,
        folderName = folderName(end-1).name;
    else
        folderName = folderName.name;
    end
    
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

    reportfig(fullfile(getenv('PROJECT'),'figures'), hfig, folderName, 'req', false,Trial.filebase,figparm,[],true,'png');
end


