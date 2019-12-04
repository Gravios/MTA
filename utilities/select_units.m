function units = select_units(varargin)
% function units = select_units(varargin)
%
% VARARGIN :
%    sessionList  -  String: {'MjgER2016'}, name use in get_session_list to call list of sessions
%    type         -  String: {'pyr','int'}, unit class selection tag
%    mode         -  String: {'get','set'}, retreive units or set model
%

% DEFARGS ----------------------------------------------------------------------
defargs = struct('sessionList',        'MjgER2016',                          ...
                 'type',               'pyr',                                ...
                 'mode',               'get'                                 ...
);
[sessionList,type,mode] = DefaultArgs(varargin,defargs,'--struct');
%-------------------------------------------------------------------------------


% MAIN -------------------------------------------------------------------------


% LOAD neuron quality data
if iscell(sessionList),
    nq = cf(@(t)  get(t.load('nq'),'nq'),  ...
            cf(@(s)  MTATrial.validate(s),  sessionList));
    nq = CatStruct(cat(1,nq{:}));    
elseif isstruct(sessionList)
    nq = cf(@(t)  get(t.load('nq'),'nq'),  ...
            af(@(s)  MTATrial.validate(s),  sessionList));
    nq = CatStruct(cat(1,nq{:}));
elseif ischar(sessionList)
    nq = cf(@(t)  get(t.load('nq'),'nq'),  ...
            af(@(s)  MTATrial.validate(s),  get_session_list(sessionList)));
    nq = CatStruct(cat(1,[nq{:}]));
elseif isa(sessionList,'MTASession'),
    sessionList.load('nq');
    nq = sessionList.nq;
else
    error('select_units:unknown input type')
end

if ~exist(fullfile(MTASession([]).path.cfg,'unit_selection_criteria.mat'),'file')&&~strcmp(mode,'set');
    % SET selection parameters
    select_units([],[],'set');        
end

switch mode
  case 'get'
% LOAD selection parameters from file    
    load(fullfile(MTASession([]).path.cfg,'unit_selection_criteria.mat'));    
    nq_type    = nq.(usp.type.fields{2})   -polyval(usp.type.pram,   nq.(usp.type.fields{1}));
    nq_quality = nq.(usp.quality.fields{2})-polyval(usp.quality.pram,nq.(usp.quality.fields{1}));
    switch type
      case 'pyr'
        units = find(nq_type<0&nq_quality>0)';
      case 'int'
        units = find(nq_type>0&nq_quality>0)';
    end

    
  case 'all'
% LOAD selection parameters from file    
    load(fullfile(MTASession([]).path.cfg,'unit_selection_criteria.mat'));    
    nq_type    = nq.(usp.type.fields{2})   -polyval(usp.type.pram,   nq.(usp.type.fields{1}));
    switch type
      case 'pyr'
        units = find(nq_type<0)';
      case 'int'
        units = find(nq_type>0)';
    end

  case 'set'
    %nqFields = {'SpkWidthR','AmpSym'};
    figH = figure(48849);clf();
    plot(nq.(nqFields{1})(nq.eDist>eDistThreshold),nq.(nqFields{2})(nq.eDist>eDistThreshold),'.')
    title('Fit for hyperplane perpendicular to selected dimensions')
    xlabel(nqFields{1});
    ylabel(nqFields{2});
    pram = draw_lines(figH,'line_fit');

    % nqFields = {'eDist','SNR'};    
    figH = figure(48849);clf();
    plot(nq.(nqFields{1})(nq.eDist>eDistThreshold),nq.(nqFields{2})(nq.eDist>eDistThreshold),'.')
    title('Fit for hyperplane perpendicular to selected dimensions')
    xlabel(nqFields{1});
    ylabel(nqFields{2});
    pramQ = draw_lines(figH,'line_fit');
    

    usp.type.pram   = pram;
    usp.type.fields = {'SpkWidthR','AmpSym'};
    usp.quality.pram   = pramQ;
    usp.quality.fields = {'eDist','SNR'};

    save(fullfile(MTASession([]).path.cfg,'unit_selection_criteria.mat'),'usp','-v7.3');
    units = [];

end
