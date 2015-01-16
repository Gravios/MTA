
% figure 1
Trial = MTASession('jg05-20120317');
xyz = Trial.load('xyz');
ang = Trial.ang.copy;
ang.create(Trial,xyz);
hfig = figure(838883);
N = histc(ang(Trial.stc{'a'},5,6,3),31:.1:37);
bar(31:.1:37,N,'histc')


vel = xyz.vel([5:9]);
vel.data(vel<0.01) = 0.01;
vel.data = log10(vel.data);

m = {'head_back','head_right'};
vbins = -.5:.05:2;
[~,mb(:,1)] = histc(vel(:,m{1}),vbins);
[~,mb(:,2)] = histc(vel(:,m{2}),vbins);

mb = MTADxyz('data',mb,'sampleRate',xyz.sampleRate);

ind = resample(Trial.stc{'t'}.cast('TimeSeries'),xyz)&nniz(mb);
%ind = nniz(mb);
%ind = ':';


A = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@mean);
B = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@std);
S = accumarray(mb(ind,:),ind(ind),repmat(numel(vbins),[1,size(mb,2)]),@sum);
A(S<100)=nan;
B(S<100)=nan;

% mean marker distance vs speed
figure,
subplot(121)
imagescnan({vbins,vbins,A'},prctile(A(nniz(A(:))),[5,95]),[],true,[0,0,0]),axis xy,
subplot(122)
imagescnan({vbins,vbins,B'},[0,2],[],true,[0,0,0]),axis xy

figure,
hist(ang(Trial.stc{'a'},'head_back','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_left',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_back','head_right',3),10:.02:37)


% subplot f


Trial = MTASession('Ed10-20140812');


xyz = Trial.load('xyz');

% lfp = Trial.lfp.copy;
% lfp.create(Trial,66);
% lfp.resample(xyz);
% 
% 
% 
% specParms = struct('nFFT',2^11,...
%                    'Fs',lfp.sampleRate,...
%                    'WinLength',2^10,...
%                    'nOverlap',2^10*.875,...
%                    'FreqRange',[1,15]);
% 
% 
% lfp.data = cat(2,lfp.data,xyz(:,1,3));
% lfp.data = cat(2,lfp.data,xyz(:,3,3));

xyz.filter(gtwin(.5,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);

figure,
plot(sq(ButFilter(xyz(:,[1:3],3)));

figure,plot(ang(:,1,3,3));


rhm = fet_rhm(Trial);
lfp.data = cat(2,lfp.data,rhm);

[ys,fs,ts,phi,fstat] = fet_spec(Trial,lfp,'mtchglong','overwrite',true);


figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy


%% Figure 5 RHM (rythmic head motion) feature versus NCP (nasal cavity pressure)
Trial = MTATrial(Ed10-20140812);

%generate features
rhm = fet_rhm(Trial,[],'wcsd');
ncp = fet_ncp(Trial,[],'wcsd');
%plot features with linked axes

figure,
sp(1) = subplot(211);
imagesc(ts,fs,log10(rhm.data)),axis xy,caxis([-6,-4])
sp(2) = subplot(212);
imagesc(ts,fs,log10(ncp.data)),axis xy,caxis([-6,-4])
linkaxes(sp,'xy');


