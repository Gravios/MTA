function processObjectPosition(fileName)
itf = c3dserver;    
openc3d(itf,0,fileName);
xyzpos = getAll3dTargets(itf);
closec3d(itf);
save([fileName '.mat'],'xyzpos');
