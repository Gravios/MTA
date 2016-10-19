function units = select_units(varargin)
% function units = select_units(varargin)
% [sessionList,eDistThreshold,type,mode,nqFields] 
% {{{'jg05-20120309','cof','all'},{'jg05-20120310','cof','all'},{'jg05-20120317','cof','all'}},...
% 15,'pyr','get',{'SpkWidthR','AmpSym'}});

% DEFARGS ----------------------------------------------------------------------
defargs = struct('sessionList',        {{'jg05-20120309.cof.all',            ...
                                         'jg05-20120310.cof.all',            ...
                                         'jg05-20120317.cof.all',            ...
                                         'Ed10-20140817.cof.all',            ...
                                         'Ed10-20140820.cof.all',            ...
                                         'ER06-20130613.cof.all',            ...
                                         'ER06-20130614.cof.all'}},          ...
                 'eDistThreshold',     15,                                   ...
                 'type',               'pyr',                                ...
                 'mode',               'get',                                ...
                 'nqFields',           {{'SpkWidthR','AmpSym'}}              ...
);%-----------------------------------------------------------------------------




[sessionList,eDistThreshold,type,mode,nqFields] = DefaultArgs(varargin,defargs,'--struct');


% MAIN -------------------------------------------------------------------------


if iscell(sessionList),
    for ses = 1:numel(sessionList),
        Trial = MTATrial.validate(sessionList{ses});
        try,
            Trial.load('nq');
        catch err
            disp(err)
            Trial.nq = NeuronQuality(Trial);
        end
        
        nq{ses} = Trial.nq;
    end
    nq = cat(1,nq{:});
    anq = CatStruct(nq);
elseif isstruct(sessionList),
    for ses = 1:numel(sessionList),
        Trial = MTATrial.validate(sessionList(ses));
        try,
            Trial.load('nq');
        catch err
            disp(err)
            Trial.nq = NeuronQuality(Trial);
        end
        
        nq{ses} = Trial.nq;
    end
    nq = cat(1,nq{:});
    anq = CatStruct(nq);
elseif ischar(sessionList)
    Trial = MTATrial(sessionList);
    Trial.load('nq');
    anq = Trial.nq;
elseif isa(sessionList,'MTASession'),
    Trial = sessionList;
    Trial.load('nq');
    anq = Trial.nq;
else
    error('select_units:unknown input type')
end


switch mode
  case 'get'
    if exist(fullfile(Trial.path.cfg,'unit_selection_criteria.mat'),'file');
        load(fullfile(Trial.path.cfg,'unit_selection_criteria.mat'));
    else
        %load(fullfile(Trial.path.cfg,'unit_selection_criteria.mat'));
    end

    nq_res = anq.(usp.fields{2})-polyval(usp.pram,anq.(usp.fields{1}));

    switch type
      case 'pyr'
        units = find(nq_res<0&anq.eDist>eDistThreshold)';
      case 'int'
        units = find(nq_res>0&anq.eDist>eDistThreshold)';
    end

  case 'set'
    figH = figure(48849);
    plot(anq.(nqFields{1})(anq.eDist>eDistThreshold),anq.(nqFields{2})(anq.eDist>eDistThreshold),'.')
    title('Fit for hyperplane perpendicular to selected dimensions')
    xlabel(nqFields{1});
    ylabel(nqFields{2});
    pram = draw_lines(figH,'line_fit');
    %saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.png');
    %saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.fig');

    usp.pram   = pram;
    usp.fields = nqFields;

    save(fullfile(Trial.path.cfg,'unit_selection_criteria.mat'),'usp','-v7.3');
    units = [];

end
