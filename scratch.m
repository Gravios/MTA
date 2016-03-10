
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,20,'low');

bfet = Trial.xyz.copy;
tsh = 15;
afet = Trial.xyz.copy;
afet.data = circshift(fxyz.data,-tsh)-circshift(fxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(fxyz(:,4,:)-fxyz(:,1,:),[1,fxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],fxyz.size(2));
% $$$ bfet.data = mean(afet(:,1:4),2)./1000.*log10(var(afet(:,1:4),[],2));
% $$$ bfet.filter('ButFilter',3,2.4,'low');
% $$$ bfet.resample(newSampleRate);

nnnfigure,plot(afet.data)
Lines(Trial.stc{'w'}(:),[],'b');

zv = afet.copy;
zv.data = log10(abs(afet(:,1))).*sign(afet(:,1));


eds = linspace(-5,5,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'w'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';


try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');


fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
man = fxyz.vel(1,[1,2]);
man.data(man.data<1e-3) =1e-3;
man.data = log10(man.data);

%edx = linspace(-0.2,1,100);
edx = linspace(-3,2,100);
edy = linspace(-5,5,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'w'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'n'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'p'};
hist2([man(ind),zv(ind)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')




ss = Trial.load('fet','3dss');
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);

figure,plot(sn)

zv = xyz.copy;
zv.data = sn;
eds = linspace(1,2,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'w'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
ind = Trial.stc{'m'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'s'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='k';


xyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));


ang = create(MTADang,Trial,xyz);
sang = [circ_dist(ang(:,1,2,1),ang(:,2,3,1)),...
        circ_dist(ang(:,2,3,1),ang(:,3,4,1)),...
        circ_dist(ang(:,3,4,1),ang(:,4,'hcom',1))];
figure,plot(abs(sum(sang,2)-circ_mean(sum(sang,2)))./ang(:,1,4,2))





av = xyz.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));
eds = linspace(0,6,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'w'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
ind = Trial.stc{'m'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'s'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='k';



fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
man = fxyz.copy;
man.data = ang(:,1,4,3).*cos(ang(:,1,4,2));

%edx = linspace(-0.2,1,100);
edx = linspace(40,160,100);
edy = linspace(0,6,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a-m'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'m'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'r'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'p'};
hist2([man(ind),zv(ind)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')




Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
xyz.
lfp = Trial.

mxyz = Trial.load('xyz');
mxyz.filter('ButFilter',3,20);
tsh = round(.15*mxyz.sampleRate);
afet = mxyz.copy;
afet.data = circshift(mxyz.data,-tsh)-circshift(mxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(mxyz(:,4,:)-mxyz(:,1,:),[1,mxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],mxyz.size(2));
afet.resample(newSampleRate);
zv = afet.copy;
zv.data = log10(abs(afet(:,1))).*sign(afet(:,1));
figure,plot(zv.data)



eds = linspace(-6,6,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'w'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
ind = Trial.stc{'m'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'s'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='k';
man = fxyz.copy;
man.data = ang(:,1,4,3).*cos(ang(:,1,4,2));

%edx = linspace(-0.2,1,100);
edx = linspace(40,160,100);
edy = linspace(0,6,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a-m'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'m'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'r'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'p'};
hist2([man(ind),zv(ind)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')




man = Trial.load('fet','lsppc');
man.filter('ButFilter',3,2,'low');


edx = linspace(-5,5,100);
edy = linspace(-5,5,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a-n'};
hist2([zv(ind,1),zv(ind,2)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'n'};
hist2([zv(ind,1),zv(ind,2)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'w'};
hist2([zv(ind,1),zv(ind,2)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'m'};
hist2([zv(ind,1),zv(ind,2)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')


fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5);
fang = create(MTADang,Trial,fxyz);
sh = 1;
sang = [circ_dist(circshift(fang(:,1,2,1),-sh),circshift(fang(:,1,2,1),sh)),...
        circ_dist(circshift(fang(:,1,3,1),-sh),circshift(fang(:,1,3,1),sh)),...
        circ_dist(circshift(fang(:,1,4,1),-sh),circshift(fang(:,1,4,1),sh)),...
        circ_dist(circshift(fang(:,1,5,1),-sh),circshift(fang(:,1,5,1),sh)),...        
        circ_dist(circshift(fang(:,1,7,1),-sh),circshift(fang(:,1,7,1),sh))...                
       ];
av = fang.copy;
av.data = sang;




figure,plot(sang)
Lines(Trial.stc{'n'}(:),[],'g');
man.data = log10(mean(abs(sang),2)./(var(sang,[],2)+1));
figure,plot(man.data)
Lines(Trial.stc{'n'}(:),[],'g');
zv = man.copy;



eds = linspace(-6,-1,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'w'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
ind = Trial.stc{'m'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'s'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='k';


tmp = [];
for s = ThreshCross(man.data,-2,1)';
    tmp(end+1) = median(man(s'));
end

tmp = [];
for s = Trial.stc{'n'}.data';
    tmp(end+1) = median(man(s'));
end

wmp = [];
for s = Trial.stc{'w'}.data';
    wmp(end+1) = median(man(s'));
end

pmp = [];
for s = Trial.stc{'p'}.data';
    pmp(end+1) = median(man(s'));
end


eds = linspace(-4,-1,100);
figure,
subplot(3,1,1),hist(tmp,eds)
subplot(3,1,2),hist(wmp,eds)
subplot(3,1,3),hist(pmp,eds)



edx = linspace(-5,-1,100);
edy = linspace(-6,6,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a-w'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'w'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'n'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'p'};
hist2([man(ind),zv(ind)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5);
fvxy = fxyz.vel([1,7],[1,2]);
fang = create(MTADang,Trial,fxyz);
figure,plot(abs(diff(fang(:,3,4,2)))),Lines(Trial.stc{'r'}(:),[], ...
                                            'r');


% FUN:ERCOR_fillgaps_RigidBody.m DATE:20160307

Session = MTASession('jg05-20120317');
xyz = Session.load('xyz');
ang = create(MTADang,Session,xyz);

bang = transform_origin(Session,xyz,'head_back','head_front',{'head_left','head_right'});
rang = transform_origin(Session,xyz,'head_left','head_right',{'head_back','head_front'});
good_index = 58320;


rb_xyz = xyz.copy;
trb_xyz = rb_xyz.copy;

rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;

nperm = factorial(rb_xyz.model.N);
bperm = 1:rb_xyz.model.N;
bperm = perms(bperm);
fperms = bperm(:,1:rb_xyz.model.N-1)';


efet = zeros([rb_xyz.size(1),nperm]);
for c = 0:NumChunks,
    if c~=NumChunks,
        ind = (c*BufferSize+1):(c+1)*BufferSize;
    else
        ind = LastPiece;
    end
    trb_xyz.data = rb_xyz(ind,:,:);
    imori = imo(trb_xyz);

    i = 1;
    for p = fperms,
        efet(ind,i) = imori(:,p(1),p(2),p(3),p(4),1);
        i = i+1;
    end
end

dtgmori = var(bsxfun(@minus,efet,efet(good_index,:)),[],2);

Efet = MTADxyz('data',efet,'sampleRate',xyz.sampleRate);

figure,imagesc(efet');
pm = 2;
figure,hist(efet(:,pm),100)
Lines(median(efet(nniz(efet),pm)),[],'r');

eds = linspace(1,1.5,1000);
figure,bar(eds,histc(efet(:,pm),eds),'histc')
Lines(median(efet(nniz(efet),pm)),[],'r');

eds = linspace(0,2.5,10000);
figure,bar(eds,histc(Efet(Session.stc{'a'},pm),eds),'histc')
Lines(median(Efet(Session.stc{'a'},pm)),[],'r');
Lines(mean(Efet(Session.stc{'a'},pm)),[],'m');
Lines(prctile(Efet(Session.stc{'a'},pm),[1,99]),[],'g');


eds = linspace(0,2.5,10000);
figure,bar(eds,histc(Efet(Session.stc{'a-e'},pm),eds),'histc')
Lines(median(Efet(Session.stc{'a-e'},pm)),[],'r');
Lines(mean(Efet(Session.stc{'a-e-s'},pm)),[],'m');
Lines(prctile(Efet(Session.stc{'a-e'},pm),[1,99]),[],'g');

figure,hist((Efet(Session.stc{'a-e'},pm)- median(Efet(Session.stc{'a-e'},pm)))./mad(Efet(Session.stc{'a-e'},pm),1),10000);

prctile(Efet(Session.stc{'a'},pm),[1,99])

bsxfun(@lt,,Efet(Session.stc{'a'},pm)


figure,hist((Efet(Session.stc{'a'},pm) - median(Efet(Session.stc{'a'},pm)))./mad(Efet(Session.stc{'a'},pm),1),10000);

mu = bsxfun(@rdivide,...
            bsxfun(@minus,...
                   Efet.data,...
                   median(Efet(Session.stc{'a'},:))),...
            mad(Efet(Session.stc{'a'},:),1));

figure,imagesc(mu');caxis([-10,10]);colormap jet
gp = fperms(:,abs(mu(125700,:))<10)


interMarDist =  imd(rb_xyz);

md = MTADxyz('data',subsref(reshape(interMarDist,size(interMarDist,1),[]),substruct('()',{':',logical(reshape(triu(ones(5,5),1),[],1))})),'sampleRate',xyz.sampleRate);
md = MTADxyz('data',bsxfun(@rdivide,...
            bsxfun(@minus,...
                   md.data,...
                   median(md(Session.stc{'a'},:))),...
            mad(md(Session.stc{'a'},:),1)),...
            'sampleRate',xyz.sampleRate);


mi = find(logical(reshape(triu(ones(5,5),1),[],1)));
mm = cell(1,2);
[mm{:}] = ind2sub([rb_xyz.model.N,rb_xyz.model.N],mi)
mm = cell2mat(mm);

gp = mm(abs(md(125700,:))>10,:)

figure,hist(md(Session.stc{'a'},1),10000);


figure,imagesc(md.data'),caxis([-10,10]);colormap jet
