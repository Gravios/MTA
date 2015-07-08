function disp_sts_sync(Trial,varargin)
if isempty(varargin)
    sts = 'rear';
else
    sts = varargin{1};
end
xyz = Trial.xyz.copy;xyz.load(Trial);
ang = Trial.ang.copy;ang.load(Trial);
vel = xyz.vel({'spine_lower','spine_middle','head_front'},[1,2]);
sp = [];
figure,
sp(1) = subplot(311);
plot(xyz(:,7,3));
Lines(Trial.stc{sts}(:),[],'r');
sp(2) = subplot(312);
plot(ang(:,3,4,3));
Lines(Trial.stc{sts}(:),[],'r');
sp(3) = subplot(313);
plot(clip(vel(:,:),0,100));
Lines(Trial.stc{sts}(:),[],'r');
linkaxes(sp,'x');