function generate_enf(path,varargin)

overwrite = false; if ~isempty(varargin) overwrite = true; end

path = '/storage/gravio/data/raw/vicon/vicon/Top Level 2/jg03/jg03-20110501/';

enfContent = [...
'[Node Information]\n'...
'TYPE=TRIAL\n'...
'Reference=\n'...
'PARENT=Session 1\n'...
'NAME=%s\n'...
'HASCHILD=FALSE\n'...
'SORTBY=1\n'...
'[TRIAL_INFO]\n'...
'CAMERACALIBRATION=\n'...
'MASK=\n'...
'[LOCK_INFO]\n'...
'CURR_USER=\n'...
];

files = dir(path)
flist = {files(~cellfun(@isempty,regexp({files.name},'.*\.Trial\d{1,5}\.enf'))).name};
clist = {files(~cellfun(@isempty,regexp({files.name},'.*\.c3d'))).name};

mtr = regexp(flist,'.*\.Trial(\d{1,5})\.enf','tokens');
for m = 1:numel(mtr), mtr{m} = mtr{m}{1}{1}; end
mtr = max(cell2mat(mtr));

for c = 1:numel(clist),    
    existingEnf = regexp(flist,[clist{c}(1:end-4),'\.Trial\d{1,5}\.enf'],...
                          'match');
    enfInd = find(~cellfun(@isempty,existingEnf),1,'first');
    if isemtpy(enfInd) enfInd = 1; end
    
    if isemtpy(existingEnf{enfInd}{1}),
        mtr = mtr+1;        
        enfFileName = [clist{c}(1:end-4) 'Trial' num2str(mtr) '.enf'];
    elseif overwrite                       
        enfFileName = existingEnf{enfInd}{1};
    else
        continue
    end
    
    fid = fopen(fullfile(path,enfFileName),'w');
    fprintf(fid,enfContent,clist{c}(1:end-4));
    fclose(fid);
end

