function feature = superSniff(Session,method,varargin)
[trialName] = DefaultArgs(varargin,{'all'});

if ~isa(Session,'MTATrial')|ischar(Session),
Trial = MTATrial(Trial,{'ang'},trialName);
end

sfet = [Trial.xyz(:,7,3),...
        Trial.ang(:,5,7,2),...
        Trial.ang(:,4,5,3),...
        Trial.ang(:,3,4,3),...
        Trial.ang(:,2,3,3)];
sfet = Filter0(gausswin(11)./sum(gausswin(11)),sfet);

num_fet = size(sfet,2);
for i = 1:num_fet
   sfet(isnan(sfet(:,i)),i) = mean(sfet(~isnan(sfet(:,i)),i));
end   

wls = WhitenSignal(sfet);
y = [];for i =1:num_fet, [y(:,:,i),f,t] = mtchglong(wls(:,i),2^8,Trial.xyzSampleRate,2^7,0.875*2^7,[],[],[],[0,20]);end


sp = [];
figure,
for i=1:num_fet,
sp(i)=subplot(num_fet,1,i);
%imagesc(t,f,log10((y(:,:,i)./repmat(max(y(:,:,i)),size(y,1),1))')),axis xy, 
%imagesc(t,f,log10(y(:,:,i)')),axis xy, 
imagesc(t,f,log10(yz(:,:,i)')),axis xy, 
end
linkaxes(sp,'xy');

figure ,imagesc(t,f,log10(sum(y,3))'),axis xy

y(y<10^-20) = 0;
figure,
for i=1:num_fet,
sp(i)=subplot(num_fet,1,i);,
hist(log10(y(y(:,10,i)~=0,10,i)),1000)
end


ym=[];
for i=1:num_fet,ym(i) = mean(mean(log10(y(y(:,1,i)~=0&~isnan(y(:,1,i)),:,i))));,end


yz = [];
for i=1:num_fet,
yz(:,:,i)=y(:,:,i).^(ym(1)-ym(i)+1);
end


figure,
%subplot(211)
imagesc(t,f,(mean(yz(:,:,1:3),3)-mean(yz(:,:,4:5),3))'),
axis xy
ys = (mean(yz(:,:,1:3),3)-mean(yz(:,:,4:5),3));
ysa = (sum(yz(:,:,1:3),3)-sum(yz(:,:,4:5),3).*1.5);
yss = yz(:,:,3)-yz(:,:,4);

Trial = Trial.load_lfp([71,75,79,83,88,92,96]);
wlfp = WhitenSignal(Trial.lfp,[],1);

yl = [];
tl = [];
fl = [];
for i =1:size(wlfp,2), 
    [yl(:,:,i),fl,tl] = mtchglong(wlfp(:,i),2^10,Trial.lfpSampleRate,2^9,0.875*2^9,[],[],[],[0,120]);
end


figure,
for i=1:4,
sp(i)=subplot(4,1,i);
if i==1,imagesc(t,f,ysa'),axis xy,caxis([0,0.0005]),
    %if i==1,imagesc(t,f,ys'),axis xy,caxis([0,0.0005]),
    %if i==1,imagesc(t,f,log10(sum(yz,3))'),axis xy,caxis([-6,-1.2]),
    %if i==1,imagesc(t,f,log10(y(:,:,1))'),axis xy,caxis([-6,-1.2]),
else  imagesc(tl,fl(fl<50),log10(yl(:,fl<50,i-1))'),axis xy,end
end
linkaxes(sp,'x');


thpow = mean(yl(:,fl>5&fl<12,3),2);
dpow = mean(yl(:,fl<5|fl>12&fl<18,3),2);
hspow = mean(yz(:,f>5&f<12),2);
lhspow = mean(yz(:,f>2&f<6),2);
hhspow = mean(yz(:,f>6&f<12),2);

[~,thi] = NearestNeighbour(tl,t,'both');

thpow = thpow(thi);
dpow = dpow(thi);



figure,
hist2([log10(thpow(hspow~=0)),log10(hspow(hspow~=0))],100,100)
figure,
hist2([log10(tvxyz(hhspow~=0&tvxyz~=0)*12),log10(hhspow(hhspow~=0&tvxyz~=0))],100,100)
figure,
hist2([log10(tvxyz(lhspow~=0&tvxyz~=0)),log10(lhspow(lhspow~=0&tvxyz~=0))],100,100)
figure,
hist2([log10(thpow(tvxyz~=0)),log10(tvxyz(tvxyz~=0))],100,100)

figure,
hist2([log10(thvxyz(hhspow~=0&thvxyz~=0)*12),log10(hhspow(hhspow~=0&thvxyz~=0))],100,100)
figure,
hist2([log10(thvxyz(lhspow~=0&thvxyz~=0)),log10(lhspow(lhspow~=0&thvxyz~=0))],100,100)
figure,
hist2([log10(thpow(thvxyz~=0)),log10(thvxyz(thvxyz~=0))],100,100)


figure,
hist2([log10(thpow(hhspow~=0)),log10(hhspow(hhspow~=0))],100,100)
figure,
hist2([log10(thpow(lhspow~=0)),log10(lhspow(lhspow~=0))],100,100)


rb = Trial.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
fxyz = Filter0(gausswin(31)./sum(gausswin(31)),Trial.com(rb));
vxyz = sqrt(sum(diff(fxyz,1).^2,2));

tx = 0:1/Trial.xyzSampleRate:size(Trial.xyz,1)/Trial.xyzSampleRate;
[~,vi] = NearestNeighbour(tx,t,'both');

tvxyz = mean([vxyz(vi),vxyz(vi+1)],2).*12;



rbh = Trial.Model.rb({'head_back','head_left','head_front','head_right'});
fhxyz = Filter0(gausswin(31)./sum(gausswin(31)),Trial.com(rbh));
vhxyz = sqrt(sum(diff(fhxyz,1).^2,2));
thvxyz = mean([vxyz(vi),vxyz(vi+1)],2).*12;


figure,
plot(log10(tvxyz(hspow~=0)),log10(hspow(hspow~=0)),'.')
figure,
hist2([log10(tvxyz(hspow~=0&tvxyz~=0)),log10(hspow(hspow~=0&tvxyz~=0))],50,50)


figure,
plot(log10(thpow(hspow~=0)),log10(tvxyz(hspow~=0)),'.')
figure,
hist2([log10(thpow(tvxyz~=0)),log10(tvxyz(tvxyz~=0))],100,100)



figure,
plot(log10(thpow(hspow~=0)./dpow(hspow~=0)),log10(tvxyz(hspow~=0)),'.')
figure,
hist2([log10(thpow(tvxyz~=0)./dpow(tvxyz~=0)),log10(tvxyz(tvxyz~=0))],100,100)


figure,
plot3(log10(tvxyz(hspow~=0)),log10(hspow(hspow~=0)),log10(thpow(hspow~=0)),'.')


wws = WhitenSignal(Trial.xyz(:,1,3));


yw = [];
[yw,fw,tw] = mtchglong(wws,2^8,Trial.xyzSampleRate,2^7,0.875*2^7,[],[],[],[0,30]);

figure,imagesc(tw,fw,log10(yw)'),axis xy


wpow = mean(yw(:,fw<6),2);
figure,
hist2([log10(wpow(tvxyz~=0&wpow~=0)),log10(tvxyz(tvxyz~=0&wpow~=0))],100,100)
figure,
hist2([log10(hhspow(tvxyz~=0&hhspow~=0)),log10(tvxyz(tvxyz~=0&hhspow~=0))],100,100)




