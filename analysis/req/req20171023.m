Trial = MTATrial.validate('jg05-20120317.cof.all');
lfp = Trial.load('lfp',65:96);

stc = Trial.load('stc','msnn_ppsvd');

xyz = preproc_xyz(Trial);

alfp = lfp.copy();
alfp.filter('RectFilter');
alfp.data = diff([RectFilter([zeros([1,size(alfp,2)]);diff(alfp.data)]);zeros([1,size(alfp,2)])]);

defspec =           struct('nFFT',2^8,'Fs',alfp.sampleRate,...
                           'WinLength',2^7,'nOverlap',2^7*.875,...
                           'FreqRange',[20,240]);

tspec =           struct('nFFT',2^11,'Fs',alfp.sampleRate,...
                         'WinLength',2^10,'nOverlap',2^10*.875,...
                         'FreqRange',[0.5,30]);


tlfp = alfp.copy();
a
ys = {};
yt = {};
for c = 12:size(lfp,2),
    tlfp.data = alfp.data(:,c);
    [ys{c},fs,ts] = fet_spec(Trial,tlfp,[],false,[],defspec);
    [yt{c},ft,tt] = fet_spec(Trial,tlfp,[],false,[],tspec);
end

figure();
subplot(311);
imagesc(ts,fs,log10(ys{5}.data'));
caxis([1.5,2.5])
% $$$ imagesc(ts,fs,nunity(ys{16}.data)');ylim([min(fs),max(fs)])
% $$$ caxis([-1,1]);
colormap('jet');
axis('xy');
subplot(312);
imagesc(tt,ft,log10(yt{5}.data'));
caxis([-2,0.5])
% $$$ imagesc(tt,ft,nunity(yt{16}.data)');ylim([min(ft),max(ft)])
% $$$ caxis([-1,2]);
colormap('jet');
axis('xy');
subplot(313);
plotSTC(Trial.stc,1);
linkaxes(findobj(gcf,'Type','axes'),'x');


yhf = cell2mat(cf(@(x) log10(x(:,22)), ys));

figure();
subplot(211);
imagesc(ts,1:32,nunity(yhf)');
colormap('jet');
caxis([-1,2]);
subplot(212);
plotSTC(Trial.stc,1);
linkaxes(findobj(gcf,'Type','axes'),'x');



yhf = cell2mat(cf(@(x) log10(x(:,30)), yt));

figure();
subplot(211);
imagesc(ts,1:32,nunity(yhf)');
colormap('jet');
caxis([-1,2]);
subplot(212);
plotSTC(Trial.stc,1);
linkaxes(findobj(gcf,'Type','axes'),'x');

yhf = ys{1}.copy();
yhf.data = nunity(cell2mat(cf(@(x) log10(x(:,25)), ys)));

yhf = yt{1}.copy();
yhf.data = nunity(cell2mat(cf(@(x) log10(x(:,15)), yt)));

rxyz = xyz.copy();
rxyz.resample(yhf);
rang = create(MTADang,Trial,rxyz);
rvxy = rxyz.vel([],[1,2]);

eds = linspace(-pi/2,pi/2,30);
bins = discretize(rang(:,'head_back','head_front',2),eds);
bins = discretize(rang(:,'spine_upper','hcom',2),eds);

bins = discretize(rang(:,'spine_upper','hcom',2),prctile(rang(nniz(rxyz),'spine_upper','hcom',2),[1:3:99]));

bins = discretize(rvxy(:,'hcom'),prctile(rvxy(nniz(rxyz),'hcom'),[1:3:99]));
sper = stc{'w+n+p+t'};
sper.resample(yhf);
sper.cast('TimeSeries');

myhfa = [];
syhfa = [];
for b = 1:32,
    myhfa(:,b) = median(yhf(logical(sper.data) & b==bins,:),'omitnan');
    syhfa(:,b) = mad(yhf(logical(sper.data) & b==bins,:),1);
end    
    
figure,imagesc(eds,1:32,myhfa)
figure,imagesc(eds,1:32,syhfa)
