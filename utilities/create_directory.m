function path = create_directory(path)
if ~exist(path,'dir'),  mkdir(path);  end