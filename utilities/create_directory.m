function path = create_directory(varargin)
path = fullfile(varargin{:});
if ~exist(path,'dir'),  mkdir(path);  end