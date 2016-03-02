function generate_enf(path,varargin)

overwrite = false; if ~isempty(varargin) overwrite = true; end

enfContent = [...
'[Node Information]\n\r'...
'TYPE=TRIAL\n\r'...
'Reference=\n\r'...
'PARENT=Session 1\n\r'...
'NAME=%s\n\r'...
'HASCHILD=FALSE\n\r'...
'SORTBY=1\n\r'...
'[TRIAL_INFO]\n\r'...
'CAMERACALIBRATION=\n\r'...
'MASK=\n\r'...
'[LOCK_INFO]\n\r'...
'CURR_USER=\n\r'...
];

files = dir(path);
flist = {files(~cellfun(@isempty,regexp({files.name},'.*\.Trial\d{1,5}\.enf'))).name};
clist = {files(~cellfun(@isempty,regexp({files.name},'.*\.c3d'))).name};

mtr = regexp(flist,'.*\.Trial(\d{1,5})\.enf','tokens');
for m = 1:numel(mtr), mtr{m} = mtr{m}{1}{1}; end
mtr = max(cellfun(@str2num,mtr));

for c = 1:numel(clist),    
    existingEnf = regexp(flist,[clist{c}(1:end-4),'\.Trial\d{1,5}\.enf'],...
                          'match');
    enfInd = find(~cellfun(@isempty,existingEnf),1,'first');
    if isempty(enfInd), enfInd = 1; existingEnf{enfInd}{1} = {};end
    
    if isempty(existingEnf{enfInd}{1}),
        mtr = mtr+1;        
        enfFileName = [clist{c}(1:end-4) '.Trial' num2str(mtr) '.enf'];
    elseif overwrite                       
        enfFileName = existingEnf{enfInd}{1};
    else
        continue
    end
    
    fid = fopen(fullfile(path,enfFileName),'w');
    fprintf(fid,enfContent,clist{c}(1:end-4));
    fclose(fid);
end

