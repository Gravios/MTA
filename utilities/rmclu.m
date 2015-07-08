function obj = rmclu(obj,clu)

cind = obj.data.clu==clu;
switch class(obj.data)
  case 'struct'
    fnames = fieldnames(obj.data);
    for f = 1:numel(fnames),
        obj.data.(fnames{f})(:,cind,:) = [];
    end
end