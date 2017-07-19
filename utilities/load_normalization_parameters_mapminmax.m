function [popMean,popStd] = load_normalization_parameters_mapminmax(functionName,varargin)

global MTA_CURRENT_PROJECT
global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
try,
    defargs = get_default_args(MTA_CURRENT_PROJECT,mfilename);
catch err
    disp(err)
    defargs = struct('referenceTrial','jg05-20120317.cof.all',...
                     'sessionList'   ,'hand_labeled',...
                     'tag'           ,'',...
                     'overwrite'     ,false);
end    
[referenceTrial,sessionList,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
if isempty(tag),
    tag = DataHash(struct('sessionList',   sessionList   ,...
                          'referenceTrial',refTrial));
end
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------
normalizationParameterFile =fullfile(MTA_PROJECT_PATH,...
                                     'analysis',...
                                     ['normalizationParameters-',functionName,'_',tag,'.mat']);

if ~exist(normalizationParameterFile,'file') || overwrite,
    Trials = af(@(Trial) MTATrial.validate(Trial)  ,get_session_list(sessionList));
    xyz    = cf(@(Trial) preproc_xyz(Trial)        ,Trials);
    fet    = cf(@(Trial,fn) feval(fn,Trial)        ,Trials, repmat({functionName},1,numel(Trials)));
    cf(@(f,t,r) f.map_to_reference_session(t,r)    ,fet,Trials,repmat({referenceTrial},1,numel(Trials)));
    for s = 1:numel(Trials), fet{s}.data(~nniz(xyz{s}),:,:,:,:) = 0;end

    % NORMALIZE feature matrix along the columns 
    zfrCat = cf(@(f) get(f,'data')    ,fet);
    zfrCat = cat(1,zfrCat{:});

    popMean = nanmean(zfrCat(nniz(zfrCat),:,:,:,:));
    popStd  = nanstd( zfrCat(nniz(zfrCat),:,:,:,:));
    save(normalizationParameterFile,'sessionList','referenceTrial','popMean','popStd');
else
    load(normalizationParameterFile);
end
% END MAIN -----------------------------------------------------------------------------------------