function out = cellstr_append_str(varargin)

if ischar(varargin{1})&&iscell(varargin{2})
    out = cellfun(@strcat,repmat({varargin{1}},size(varargin{2})),varargin{2},'UniformOutput',false);
elseif ischar(varargin{2})&&iscell(varargin{1})
    out = cellfun(@strcat,varargin{1},repmat({varargin{2}},size(varargin{1})),'UniformOutput',false);
else
    for a = 1:numel(varargin),
        disp(['arg',num2str(a),': ',class(varargin{a})])
    end
    error('MTA:utilities:cellstr_append_str, Unknow input pattern')
end