Trial = MTATrial.validate('Ed05-20140529.ont.all');
stc = Trial.load('stc','hand_labeled_rev1_Ed');
Trial = MTATrial.validate('jg05-20120317.cof.all');
stc = Trial.load('stc','hand_labeled_rev3_jg');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
hang = Trial.transform_origin;



edx = linspace(35,57,100);
ind = stc{'n'};
figure,
bar(edx,histc(ang(ind,3,4,3),edx),'histc');

edy = linspace(-1,0.2,100);
figure,
bar(edy,histc(ang(ind,3,4,2),edy),'histc');

edy = linspace(0,pi/2,100);
figure,
bar(edy,histc(ang(ind,1,2,2),edy),'histc');




figure,
hist2(sq(ang(ind,3,4,2:3)),edy,edx)



edx = linspace(0.2,.9,100);
edy = linspace(60,95,100);
figure,
hist2([ang(ind,1,3,2),xyz(ind,4,3)],edx,edy)


ind = stc{'n'}(:,1);
ind = [ind-20,ind+60];
ind(1)=[];
edx = linspace(0.2,.9,100);
edy = linspace(60,95,100);
figure,
hist2([ang(ind,1,3,2),xyz(ind,4,3)],edx,edy)