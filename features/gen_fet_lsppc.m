function gen_fet_lsppc(Trial);

xyz = Trial.load('xyz');
afet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-5)-circshift(xyz(:,:,[1,2]),5);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
afet.data = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));



mag = zeros([afet.size(1),1]);

% GET avalible markers from list
mrkInd = Trial.xyz.model.gmi({'spine_lower','pelvis_root','spine_middle','spine_upper',...
                              'head_back', 'head_front'});
mrkInd(~mrkInd) = [];
% REMOVE unavailable markers

for i= 1:afet.size(1),
    mag(i) = PPC(afet(i,mrkInd));
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
