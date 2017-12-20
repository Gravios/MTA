function Session = load(Session,varargin)
%Session = load(Session,varargin)
%load the session file
%
%  Session:     MTASession, The data object which synchronizes and
%                           holds all experimental data.
%
%  field:       string,     The name of a field which belongs to
%                           session, which will be loaded from the
%                           data file and synchronized with the
%                           current sync.
%
fvarargin = {};
if numel(varargin)>1,
    fvarargin = varargin(2:end);
end
[field] = DefaultArgs(varargin,{[]});
if ~isempty(field),
    pattern = {'MTAData';'MTASpk';'MTAStateCollection'};
    F_classes = cat(2,superclasses(Session.(field))',{class(Session.(field))});
    if any(subsref(~cellfun(@isempty,regexp(repmat(F_classes,[numel(pattern),1]),...
                                            repmat(pattern,[1,numel(F_classes)]))),...
                   substruct('()',{':'}))),
        %if isa(Session.(field),'MTAData')||isa(Session.(field),'MTASpk'),      
        if nargout==1,
            Data = Session.(field).copy;
            Data = Data.load(Session,fvarargin{:}); 
            Session = Data;
        else
            Session.(field).load(Session,fvarargin{:}); 
        end
    else
        switch varargin{1}
          case 'nq'
            ds = load(fullfile(Session.spath, [Session.name '.NeuronQuality.mat']));
            Session.nq = ds.nq;
        end
    end
else
    load(fullfile(Session.spath, [Session.filebase '.ses.mat']));    
end%if emptyfield

