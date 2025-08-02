function path = generate_mta_paths(meta)
% function path = generate_mta_paths(meta)
% create all data directiories and store their paths
MTA_DATA_PATH = getenv('MTA_DATA_PATH');
path.raw.data        = create_directory(fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase));
path.raw.ephys       = create_directory(fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase,meta.sessionName));
path.processed.ephys = create_directory(fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName));

create_directory(fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName));
create_directory(fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName));

for maze = meta.mazeList
    maze = maze{1};
    path.raw.xyz.(maze)       = create_directory(fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName,maze));
    path.processed.xyz.(maze) = create_directory(fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName,maze));
end