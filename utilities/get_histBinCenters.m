function [varargout] = get_histBinCenters(varargin)

if iscell(varargin{1}),
    edgs = varargin{1};
else    
    edgs = varargin;
end

edgs = cellfun(@plus,edgs,...
            cellfun(@rdivide,...
               cellfun(@diff,...
                  cellfun(@subsref,edgs,repmat({substruct('()',{[1,2]})},size(edgs)),'uniformoutput',false),'uniformoutput',false),...
            repmat({2},size(edgs)),'uniformoutput',false),'uniformoutput',false);

% $$$ edgs = cellfun(@subsasgn,edgs,cellfun(@substruct,repmat({'()'},size(edgs)),cellfun(@mat2cell,cellfun(@length,edgs,'uniformoutput',false),repmat({1},size(edgs)),repmat({1},size(edgs)),'uniformoutput',false),'uniformoutput',false),{[],[]},'uniformoutput',false);
for e = 1:numel(edgs),
    edgs{e}(end) = [];
end

varargout = edgs;
if nargout<numel(varargout),
varargout = varargout(1:nargout);
end