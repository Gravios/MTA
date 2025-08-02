% req20211210

% CSD type stuff gamma frequency

configure_default_args();
MjgER2016_load_data();

Trial = Trials{1};
sampleRate = 30;

xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
stateColors = 'krggbbmy';


lfp = Trial.load('lfp',[57,60,61,64]);
lfp = Trial.load('lfp',[57,58,60,61,63,64]);
lfp = Trial.load('lfp',[57:64]);
lfp = Trial.load('lfp',[1,4,5,8]);
lfp = Trial.load('lfp',[1,2,4,5,7,8]);
lfp.filter('ButFilter',4,[30],'low');
lfp.filter('ButFilter',4,[300],'low');
% $$$ 
% $$$ figure,
% $$$ plot(nunity(lfp(:,1)))
% $$$ hold('on');
% $$$ plot(nunity(diff(lfp(:,:),1,2)));

clfp = lfp.copy();
clfp.data = [nunity(lfp(:,:)),nunity(diff(lfp(:,:),1,2))];

specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  clfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[lys,lfs,lts] = fet_spec(Trial,clfp,[],false,[],specArgsTheta);




figure();
sax = gobjects([0,1]);
sax(end+1) = subplot(411);
hold('on');
imagesc(lts,lfs,log10(lys(:,:,1))'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(412);
hold('on');
imagesc(lts,lfs,log10(lys(:,:,2))'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(413);
hold('on');
imagesc(lts,lfs,log10(lys(:,:,3))'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(414);
hold(sax(end),'on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
xlabel('Time (s)');
linkaxes(sax,'x');


clfp = lfp.copy();
clfp.data = [nunity(lfp(:,:)),nunity(diff(lfp(:,:),1,2))];
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  clfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[lys,lfs,lts] = fet_spec(Trial,clfp,[],false,[],specArgsTheta);

clfp = lfp.copy();
clfp.data = [nunity(lfp(1:1e5,:)),nunity(diff(lfp(1:1e5,:),1,2)),nunity(diff(lfp(1:1e5,[1,end]),1,2))];
clfp.data = [nunity(diff(lfp(1:2e5,:))),nunity(diff(diff(lfp(1:2e5,:)),1,2)),nunity(diff(diff(lfp(1:2e5, ...
                                                  [1,end])),1,2))];
ind = 2e5;
clfp.data = [(diff(lfp(1:ind,[1,4,end]))),...
             (diff(lfp(1:ind,1))-diff(lfp(1:ind,2))),...
             (diff(lfp(1:ind,1))-diff(lfp(1:ind,3))),...             
             (diff(lfp(1:ind,1))-diff(lfp(1:ind,4))),...                          
             (diff(lfp(1:ind,2))-diff(lfp(1:ind,4))),...                                       
             (diff(lfp(1:ind,3))-diff(lfp(1:ind,4))),...                                                    
             (diff(lfp(1:ind,3))-diff(lfp(1:ind,5))),...                                                    
             (diff(lfp(1:ind,4))-diff(lfp(1:ind,6))),...                                                    
             (diff(diff(lfp(1:ind,[1,end])),1,2))];
specArgsTheta = struct('nFFT',2^8,...
                  'Fs',  clfp.sampleRate,...
                  'WinLength',2^7,...
                  'nOverlap',2^7*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[25,180]);
[gys,gfs,gts] = fet_spec(Trial,clfp,[],false,[],specArgsTheta);



figure();
sax = gobjects([0,1]);
sax(end+1) = subplot(711);
hold('on');
imagesc(gts,gfs,RectFilter(RectFilter(log10(gys(:,:,8))',3,3)',3,3)'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(712);
hold('on');
imagesc(gts,gfs,RectFilter(RectFilter(log10(gys(:,:,4))',3,3)',3,3)'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(713);
hold('on');
imagesc(gts,gfs,RectFilter(RectFilter(log10(gys(:,:,5))',3,3)',3,3)'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(714);
hold('on');
imagesc(gts,gfs,RectFilter(RectFilter(log10(gys(:,:,7))',3,3)',3,3)'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(715);
hold('on');
imagesc(gts,gfs,RectFilter(RectFilter(log10(gys(:,:,9))',3,3)',3,3)'); 
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(716);
hold('on');
imagesc(gts,gfs,nunity(RectFilter(RectFilter(log10(gys(:,:,10))',3,3)',3,3))'); 
plot((1:size(lfp,1))./lfp.sampleRate,nunity(lfp(:,1)-lfp(:,4))*20+60);
plot((1:size(lfp,1))./lfp.sampleRate,nunity(lfp(:,1))*20+60);
axis('tight'); axis('xy'); colormap('jet');
sax(end+1) = subplot(717);
hold(sax(end),'on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
xlabel('Time (s)');
linkaxes(sax,'x');

msp = [];
for c = [5,7,9,10]
msp(:,:,end+1) = RectFilter(RectFilter(log10(gys(:,:,[c]))',3,3)',3,3)';
end

figure,imagesc(gts,gfs,1./(mean(msp,3)./std(msp,[],3))); 
axis('tight'); axis('xy'); colormap('jet');
% Compute new theta states 
% Compare and rewiew labeling reliabilty with spatial gradient measure
frin = [6,12];
frout = reshape([1,4,12,14],2,[])';

thfin = logical(WithinRanges(s.f,frin));
thfout = logical(WithinRanges(s.f,frout));

newPeriods = HmmStateSegment(gCheckEegStates.FileBase, detPer,frin,frout,ch2use)/gCheckEegStates.eFs;