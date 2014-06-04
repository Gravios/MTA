function units = select_units(varargin)
% $$$ function units = select_units(varargin)
% $$$ [SessionList,eDistThreshold,type,mode,nqFields] 
% $$$ {{{'jg05-20120309','cof','all'},{'jg05-20120310','cof','all'},{'jg05-20120317','cof','all'}},...
% $$$ 15,'pyr','get',{'SpkWidthR','AmpSym'}});

[SessionList,eDistThreshold,type,mode,nqFields] = DefaultArgs(varargin,...
{{{'jg05-20120309','cof','all'},{'jg05-20120310','cof','all'},{'jg05-20120317','cof','all'}},...
15,'pyr','get',{'SpkWidthR','AmpSym'}});

if iscell(SessionList),
    for ses = 1:numel(SessionList),
        Trial = MTATrial(SessionList{ses}{1},SessionList{ses}{3},SessionList{ses}{2});
        Trial.load('nq');
        nq{ses} = Trial.nq;
    end
    nq = cat(1,nq{:});
    anq = CatStruct(nq);
elseif ischar(SessionList)
    Trial = MTATrial(SessionList);
    Trial.load('nq');
    anq = Trial.nq;

elseif isa(SessionList,'MTASession'),
    Trial = SessionList;
    Trial.load('nq');
    anq = Trial.nq;
else
    error('select_units:unknown input type')
end

switch mode
  case 'get'
    load(fullfile(Trial.path.MTAPath,'unit_selection_criteria.mat'));
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

    save(fullfile(Trial.path.MTAPath,'unit_selection_criteria.mat'),'usp','-v7.3');
    units = [];

end
