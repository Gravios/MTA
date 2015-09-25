function gen_fet_lsppc(Trial);

xyz = Trial.load('xyz');
afet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-5)-circshift(xyz(:,:,[1,2]),5);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
afet.data = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));



mag = zeros([afet.size(1),1]);
for i= 1:afet.size(1),
mag(i) = PPC(afet(i,[1:5,7]));
end

msync = Trial.xyz.sync.copy;
msync.data = msync.sync.data;


man = MTADfet(Trial.spath,Trial.filebase,...
               mag,...
               Trial.xyz.sampleRate,...
               msync,...
               msync(1),...
               [],[],[],'lower_spine_trajectory_yaw_PPC','lsppc','y');
man.save;
