function [popPS] = load_normalization_parameters_mapminmax(functionName,varargin)

global MTA_CURRENT_PROJECT
global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
try,
    defargs = get_default_args(MTA_CURRENT_PROJECT,mfilename);
catch err
    disp(err);
    for s = 1:numel(err.stack)
        disp(err.stack(s));
    end
    
    defargs = struct('referenceTrial','jg05-20120317.cof.all',                                   ...
                     'sessionList'   ,'hand_labeled',                                            ...
                     'normalize'     ,true,                                                      ...
                     'tag'           ,'',                                                        ...
                     'overwrite'     ,false);
end    
[referenceTrial,sessionList,normalize,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
if isempty(tag),
    tag = DataHash(struct('sessionList',   sessionList   ,...
                          'referenceTrial',referenceTrial,...
                          'normalize',     normalize));
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
normalizationParameterFile =fullfile(MTA_PROJECT_PATH,...
                                     'analysis',...
                                     ['normalizationParameters-mapminmax-',functionName,'_',tag,'.mat']);

if ~exist(normalizationParameterFile,'file') || overwrite,
    Trials = af(@(Trial) MTATrial.validate(Trial)  ,get_session_list(sessionList));
    xyz    = cf(@(Trial) preproc_xyz(Trial)        ,Trials);
    fet    = cf(@(Trial,fn) feval(fn,Trial)        ,Trials, repmat({functionName},1,numel(Trials)));
    cf(@(f,t,r) f.map_to_reference_session(t,r)    ,fet,Trials,repmat({referenceTrial},1,numel(Trials)));
    for s = 1:numel(Trials), fet{s}.data(~nniz(xyz{s}),:,:,:,:) = 0;end

    if normalize,
        [refMean,refStd] = load_normalization_parameters_unity(functionName,...
                                                               referenceTrial,...
                                                               sessionList,...
                                                               [],...
                                                               overwrite);
        cf(@(f,m,s) f.unity(@nan,m,s),...
            fet,...
            repmat({refMean},1,numel(Trials)),...
            repmat({refStd}, 1,numel(Trials)));
    end
    

% $$$     [StcRnd,labelingEpochs,trainingFeatures] = cf(@(s,f,sts) ...
% $$$             resample_whole_state_bootstrap_noisy_trim(s,f,sts,[100],10000),...
% $$$             Stc,fet,repmat({states},1,numel(Trials)));
    
    mm = cf(@(f)  bsxfun(@plus,prctile(f.data,[1,99]),[-4;4])    ,fet);
    mm = cat(3,mm{:});    
    mm = [min(mm(1,:,:),[],3)',max(mm(2,:,:),[],3)'];
% $$$     fdata = [];
% $$$     for s = 1:numSessions,
% $$$         fdata = cat(1,fdata,trainingFeatures{s}.data);
% $$$     end
% $$$     mm = bsxfun(@plus,prctile(fdata,[0.1,99.9]),[-2;2]);
  
    popPS.name   = 'mapminmax';
    popPS.xrows  = size(fet{1},2);
    popPS.xmax   = mm(:,2);
    popPS.xmin   = mm(:,1);
    popPS.xrange = diff(mm,1,2);
    popPS.yrows  = size(fet{1},2);
    popPS.ymax   = 1;
    popPS.ymin   =-1;
    popPS.yrange = 2;
    popPS.no_change = 0;
    popPS.gain      = 2./popPS.xrange;
    popPS.xoffset   = popPS.xmin;

    % NORMALIZE feature matrix along the columns 
    save(normalizationParameterFile,'sessionList','referenceTrial','normalize','popPS');
else
    load(normalizationParameterFile);
end
% END MAIN -----------------------------------------------------------------------------------------