

Trial = MTATrial('jg05-20120317');

% XYZ Positions of Markers
xyz = Trial.load('xyz');

% COM Body Center of Mass
rbb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},xyz.com(rbb));




% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

% COM Body Center of Mass
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

fbcom = bcom.copy;
fbcom.filter('ButFilter',3,2.5,'low');

% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% FVEL Filtered marker speeds in XY plane
fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0)=.1;
fvel.data = log10(fvel.data);


llist = {};
figure,hold on

% $$$ llist{end+1} = 'fbcom';               % Filtered Center of mass of body 
% $$$ plot(nunity(fbcom(:,3)))
% $$$ 
% $$$ llist{end+1} = 'body pitch';          % Pitch of marker set {'spine_lower','spine_upper'} 
% $$$ plot(nunity(ang(:,1,4,2)))    

llist{end+1} = 'upper spine pitch';   % Pitch of marker set {'spine_middle','spine_upper'} 
plot(nunity(fang(:,3,4,2)))
llist{end+1} = 'd(upper spine pitch)/dt';   % Pitch of marker set {'spine_middle','spine_upper'} 
plot(nunity(diff(fang(:,3,4,2))))

llist{end+1} = 'head height';         % displacement in z-axis relative to floor
plot(nunity(fxyz(:,5,3)))
llist{end+1} = 'd(head height)/dt'    % first derivative of low passed filtered head hight
plot(nunity(diff(fxyz(:,5,3))))



% $$$ llist{end+1} = 'd(fbcom)/dt'    % first derivative of low passed filtered head hight
% $$$ plot(nunity(diff(fbcom(:,3))))

legend(llist)


figure,plot(cos(fang(:,3,4,2)).*fang(:,3,4,3))
Lines(Trial.stc{'r'}(:),[],'r');

figure,hold on,m2 = 4;
eds = linspace(20,250,100);
eds = linspace(20,180,100);
ind = Trial.stc{'a-r-m'};
hn = bar(eds,histc(cos(fang(ind,1,m2,2)).*fang(ind,1,m2,3),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha =  0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(cos(fang(ind,1,m2,2)).*fang(ind,1,m2,3),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha =  0;
ind = Trial.stc{'m'};
hs = bar(eds,histc(cos(fang(ind,1,m2,2)).*fang(ind,1,m2,3),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;
hs.EdgeAlpha =  0;


figure,
eds = linspace(20,170,100);
ads = linspace(-0.8,1.5,100);
ind = Trial.stc{'a'};
hist2([fang(ind,3,4,2),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,500])

figure,% Body length xy VS lower spine marker speed
eds = linspace(20,170,100);
ads = linspace(-0.5,2,100);
ind = Trial.stc{'n'};
hist2([fvel(ind,1),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,200])

figure,
eds = linspace(20,170,80);
ads = linspace(-5,.5,80);
ind = Trial.stc{'m'};
hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-20)))),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,200])



vangp = MTADxyz('data',log10(abs(nunity(diff(fang(:,3,4,2))))),'sampleRate',xyz.sampleRate);

figure,
ads = linspace(-4,1,100);
eds = linspace(20,170,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])


figure,
ads = linspace(-4,1,100);
eds = linspace(-1,1.8,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])

figure,
ads = linspace(-4,1,100);
eds = linspace(60,260,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])


figure,
ads = linspace(20,170,100);
eds = linspace(60,260,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])






figure,hold on
eds = linspace(-.5,2,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(ang(ind,1,4,2),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
ind = Trial.stc{'r'};
hs = bar(eds,histc(ang(ind,1,4,2),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;




%% Detect Rearing Periods

figure,plot(nunity(bcom(:,3)))


v = xyz.vel(1:8,[1,2]);


%% Find best low pass filter for detecting walk periods
% $$$ dp = []
% $$$ frange = .5:.2:4;
% $$$ for i = frange,
% $$$ vt = v.copy;
% $$$ vt.filter('ButFilter',3,i,'low');
% $$$ vt.data(vt.data<0) = .001;
% $$$ vta = log10(vt(Trial.stc{'a-w'},:));
% $$$ vtw = log10(vt(Trial.stc{'a&w'},:));
% $$$ nw = nniz(vtw);
% $$$ na = nniz(vta);
% $$$ dp(end+1,:) = (mean(vtw(nw,:))-mean(vta(na,:)))./(.5*sqrt(var(vtw(nw,:))+var(vta(na,:))));
% $$$ end


% $$$ [~,fid] = max(dp(:,1));

% $$$ v.filter('ButFilter',3,frange(fid),'low');
v.filter('ButFilter',3,2.5,'low');
% $$$ figure,plot(v.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');

v.data(v.data<.1) = .1;
v.data = log10(v.data);

pv =MTADxyz('data',prod(v(:,[1,4,7]),2),'sampleRate',xyz.sampleRate);
pv.data(pv<0)= 1e-15;
pv.data = log10(pv.data);
figure,plot(pv.data);
Lines(Trial.stc{'w'}(:),[],'b');



figure,hold on
eds = linspace(-15,2.5,100);
ind = Trial.stc{'a-w-r-n'};
bar(eds,histc(pv(ind,1),eds),'histc');
ind = Trial.stc{'w+n'};
h = bar(eds,histc(pv(ind,1),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;


%lb pitch vs lb vel
peds = linspace(-10,2.5,100);
veds = linspace(-.5,2,100);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([ang(ind,1,2,2),v(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w+n'};
figure,hist2([ang(ind,1,2,2),v(ind,1)],peds,veds);caxis([0,500])

% 
peds = linspace(110,170,100);
veds = linspace(-10,2.5,100);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([ang(ind,1,4,3),pv(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w+n'};
figure,hist2([ang(ind,1,4,3),pv(ind,1)],peds,veds);caxis([0,500])


peds = linspace(100,160,70);
veds = linspace(-2,1,70);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'n'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,200])



peds = linspace(-4,0,70);
ind = Trial.stc{'n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'w'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'a-w-r-n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])

veds = linspace(100,160,70);
peds = linspace(-4,0,70);
ind = Trial.stc{'n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,2,1),circshift(fang(ind,1,2,1),-20)))),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])
ind = Trial.stc{'w'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,2,1),circshift(fang(ind,1,2,1),-20)))),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])
ind = Trial.stc{'a-w-r-n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,2,1),circshift(fang(ind,1,2,1),-20)))),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])



figure,hold on
eds = linspace(-.5,2,100);
ind = Trial.stc{'a-w-r'};
bar(eds,histc(v(ind,1),eds),'histc');
ind = Trial.stc{'w'};
h = bar(eds,histc(v(ind,1),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;

figure,plot(v(:,1).^2./abs(v(:,1)-v(:,7)));
Lines(Trial.stc{'w'}(:),[],'m');

vco = v.data;
vco = bsxfun(@minus,v.data,mean(v.data)


v.data(v.data<0) = .001;
v.data = log10(v.data);






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


xyz.filter('ButFilter',3,2,'low');
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





ind = Trial.stc{'m'};
figure,
hist2([a(ind,1,4,3),a(ind,2,7,3)],linspace(110,170,70),linspace(70,240,70));
colormap jet
caxis([0,300]);

caxis([0,800]);

ind = Trial.stc{'a-m'};
figure,
hist2([a(ind,1,4,3),a(ind,4,7,3)],linspace(100,170,70),linspace(70,170,70));
colormap jet
caxis([0,300]);

ind = Trial.stc{'m'};
figure,
hist2([a(ind,1,4,3),circ_dist(a(ind,1,2,1),a(ind,1,4,1))],linspace(100,170,70),linspace(-pi/2,pi/2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);



figure,hold on
eds = linspace(100,160,100);
ind = Trial.stc{'a-m-s'};
bar(eds,histc(a(ind,1,4,3),eds),'histc');
ind = Trial.stc{'m'};
h = bar(eds,histc(a(ind,1,4,3),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;




figure,hold on
eds = linspace(-pi,pi,100);
ind = Trial.stc{'a-m'};
bar(eds,histc(circ_dist(a(ind,1,2,1),a(ind,1,4,1)),eds),'histc');
ind = Trial.stc{'m'};
h = bar(eds,histc(circ_dist(a(ind,1,2,1),a(ind,1,4,1)),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;



