

Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');

v = xyz.vel(1:8,[1,2]);

%% Find best low pass filter for detecting walk periods
dp = []
frange = .5:.2:4;
for i = frange,
vt = v.copy;
vt.filter('ButFilter',3,i,'low');
vt.data(vt.data<0) = .001;
vta = log10(vt(Trial.stc{'a-w'},:));
vtw = log10(vt(Trial.stc{'a&w'},:));
nw = nniz(vtw);
na = nniz(vta);
dp(end+1,:) = (mean(vtw(nw,:))-mean(vta(na,:)))./(.5*sqrt(var(vtw(nw,:))+var(vta(na,:))));
end


[~,fid] = max(dp(:,1));

v.filter('ButFilter',3,frange(fid),'low');
% $$$ figure,plot(v.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');

v.data(v.data<.1) = .1;
v.data = log10(v.data);



figure,hold on
eds = linspace(-.5,2,100);
ind = Trial.stc{'a-w-r-n'};
bar(eds,histc(v(ind,1),eds),'histc');
ind = Trial.stc{'n'};
h = bar(eds,histc(v(ind,1),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;

figure,plot(v(:,1).^2./abs(v(:,1)-v(:,7)));
Lines(Trial.stc{'w'}(:),[],'m');

vco = v.data;
vco = bsxfun(@minus,v.data,mean(v.data)


v.data(v.data<0) = .001;
v.data = log10(v.data);





rbb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
bcom = MTADfet(Trial.spath,...
               [],...
               sq(xyz.com(rbb)),...
               xyz.sampleRate,...
               xyz.sync.copy,...
               xyz.origin,...
               [],[],[],...
               'bodyCOM',...
               'bcom',...
               'b');

[bys,bfs,bts,bphi,bfst] = fet_spec(Trial,bcom);

figure,imagesc(bts,bfs,log10(bys(:,:,1,1))'),axis xy,colormap jet

figure,
for i = 1:bys.size(3),
sp(i) = subplot(bys.size(3),1,i);
imagesc(bts,bfs,log10(bys(:,:,i,i))')
axis xy
colormap jet
end
linkaxes(sp,'xy');
Lines(Trial.stc{'w',1}(:),[],'b');
Lines(Trial.stc{'r',1}(:),[],'r');


figure,imagesc(bts,bfs,bys(:,:,1,3)'),axis xy, colormap jet
figure,imagesc(bts,bfs,bys(:,:,2,3)'),axis xy, colormap jet

figure,imagesc(bts,bfs,(bys(:,:,1,3).*bys(:,:,1,2).*bys(:,:,2,3))'),axis xy, colormap jet


figure,imagesc(bts,bfs,(bys(:,:,1,3)+bys(:,:,1,2)+bys(:,:,2,3))'),axis xy, colormap jet
caxis([2,3])
Lines(Trial.stc{'m',1}(:),[],'m');
Lines(Trial.stc{'k',1}(:),[],'r');

figure,imagesc(bts,bfs,log10(bfst(:,:,3))'),axis xy, colormap jet

figure,plot(bcom(:,1,3))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'r'}(:),[],'r');


xyz.filter('ButFilter',3,3,'low');
a = create(MTADang,Trial,xyz);
h = Trial.transformOrigin(xyz,'head_back','head_front',{'head_left','head_right'});

r = xyz.copy;
r.data = h.roll;


%Groom
ind = Trial.stc{'m'};
figure,
hist2([circ_dist(a(ind,1,3,1),a(ind,4,7,1)),r(ind)-.2],linspace(-pi,pi,70),linspace(-pi/2,pi/2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);

ind = Trial.stc{'w'};
figure,
hist2([a(ind,1,4,3),v(ind,1)],linspace(110,170,70),linspace(-.5,2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);





ind = Trial.stc{'a-m'};
figure,
hist2([a(ind,1,4,3),a(ind,1,7,3)],linspace(110,170,70),linspace(70,280,70));
colormap jet
caxis([0,300]);

caxis([0,800]);




