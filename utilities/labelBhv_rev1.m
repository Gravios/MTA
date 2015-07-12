

Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');

v = xyz.vel(1:8,[1,2]);

dp = []
for i = .5:.2:4
vt = v.copy;
vt.filter('ButFilter',3,i,'low');
vt.data(vt.data<0) = .001;
vta = log10(vt(Trial.stc{'a-w'},:));
vtw = log10(vt(Trial.stc{'a&w'},:));
dp(end+1,:) = (mean(vta)-mean(vtw))./(.5*sqrt(var(vta)+var(vtw)));
end


v.filter('ButFilter',3,2,'low');
figure,plot(v.data)
Lines(Trial.stc{'w'}(:),[],'b');

figure,plot(v(:,1).^2./abs(v(:,1)-v(:,7)));
Lines(Trial.stc{'w'}(:),[],'m');

vco = v.data;
vco = bsxfun(@minus,v.data,mean(v.data)


v.data(v.data<0) = .001;
v.data = log10(v.data);



