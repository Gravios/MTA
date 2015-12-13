function out = save(Data,varargin)        
out = false;
warning('off','MATLAB:structOnObject');
Dstruct = struct(Data);
warning('on','MATLAB:structOnObject');
if isempty(varargin),
    sargs = {Data.fpath,'-struct','Dstruct','-v7.3'};
else
    sargs = {varargin{1},'-struct','Dstruct','-v7.3'};
end
save(sargs{:});
out = true;

end
