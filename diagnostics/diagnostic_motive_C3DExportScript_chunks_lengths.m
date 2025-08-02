itf = c3dserver;
fileloc = 'E:\motive\rs06-20151201-sof\vicon';
files = dir(fileloc);
c3dlist = {files(~cellfun(@isempty,regexp({files.name},'Trial001_\d*_vicon.c3d'))).name};
       
xs = zeros(size(c3dlist));
for c = 1:numel(c3dlist)
openc3d(itf,0,fullfile(fileloc,c3dlist{c}));
[xyz,~] = getAll3dTargets(itf);
xs(c) = size(xyz,1);
closec3d(itf)
end

itf.release


