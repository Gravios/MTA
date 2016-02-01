
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

figure,plot(afet.data)
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

