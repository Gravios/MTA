% req20201111
%    Tags: interneuron ratemaps
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: BhvPlaceCode
%    Description: spectral decomposition of marker kinematics
%                 


Trial = MTATrial.validate('jg05-20120317.cof.all');
fet = fet_bref(Trial);
ofet = fet_bref(Trial);
stc = Trial.load('stc','hand_labeled_rev3_jg');

figure();
hold('on');
plot(ofet(:,30))
Lines(stc{'r'}.data(:,1),[],'r');
Lines(stc{'r'}.data(:,2),[],'k');

fet.data = ofet(:,[17:4:25,26:2:30]);
fet.data = ofet(:,[16:end]);
fet.data = ofet(:,[16,17,20,21,24,25,26,28,30]);
fet.data = ofet(:,[17]);
defspec = struct('nFFT',2^8,'Fs',fet.sampleRate,...
                 'WinLength',2^8,'nOverlap',2^8-2^4,...
                 'FreqRange',[1,20]);

% $$$ defspec = struct('nFFT',2^8,'Fs',fet.sampleRate,...
% $$$                  'WinLength',2^7,'nOverlap',2^7-2^4,...
% $$$                  'FreqRange',[1,20]);


[ys,fs,ts] = fet_spec(Trial,fet,'defspec',defspec,'flagCrossSpec',false);

mYs = median(ys(stc{'a'},:,:));
sYs = std(ys(stc{'a'},:,:));
uYs = ys.copy();
uYs.data = nunity(ys.data,@nan,mYs,sYs);

mode = 'log10';
hfig = figure();
sp = gobjects([0,1]);
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hfig.Position(3)    = 29.7;
hfig.Position(4)    = 21.0;
sp(end+1) = subplot(411);
plot((1:size(fet.data,1))./fet.sampleRate,fet(:,:));
sp(end+1) = subplot(412);
c = 1;
switch mode
  case 'unity', imagesc(ts,fs,uYs(:,:,c,1)');
  case 'funity', imagesc(ts,fs,bsxfun(@rdivide,uYs(:,:,c,1)+0.5,sum(uYs(:,:,c,1)+0.5,2))');    
  case 'log10', imagesc(ts,fs,log10(ys(:,:,c,1))');
  case 'fnorm', imagesc(ts,fs,bsxfun(@rdivide,ys(:,:,c,1),sum(ys(:,:,c,1),2))');
end
axis('xy');
sp(end+1) = subplot(413);
c = 3;
switch mode
  case 'unity', imagesc(ts,fs,uYs(:,:,c,1)');
  case 'funity', imagesc(ts,fs,bsxfun(@rdivide,uYs(:,:,1,1)+0.5,sum(uYs(:,:,c,1)+0.5,2))');    
  case 'log10', imagesc(ts,fs,log10(ys(:,:,c,1))'); 
  case 'fnorm', imagesc(ts,fs,bsxfun(@rdivide,ys(:,:,1,1),sum(ys(:,:,c,1),2))');
end
axis('xy');
sp(end+1) = subplot(414);
c = 4;
switch mode
  case 'unity', imagesc(ts,fs,uYs(:,:,c,1)');
  case 'funity', imagesc(ts,fs,bsxfun(@rdivide,uYs(:,:,1,1)+0.5,sum(uYs(:,:,c,1)+0.5,2))');    
  case 'log10', imagesc(ts,fs,log10(ys(:,:,c,1))');
  case 'fnorm', imagesc(ts,fs,bsxfun(@rdivide,ys(:,:,1,1),sum(ys(:,:,c,1),2))');
end
axis('xy');
linkaxes(sp,'x');
colormap('jet');
switch mode,
  case 'unity' , ForAllSubplots('caxis([-0.2,5]);');
  case 'funity', ForAllSubplots('caxis([0.001,0.1]);');
  case 'fnorm' , ForAllSubplots('caxis([0,0.1]);');
  case 'log10' , ForAllSubplots('caxis([-4,0]);');
end


stcm = ys.copy();
stcm.data = stc2mat(stc,ys,{'rear','walk','turn','pause','groom','sit'});

ffet = ofet.copy();
ffet.filter('ButFilter',4,4,'low');
ffet.resample(ys);

fys = ys.copy();
fys.data = reshape(fys(:,1:2:end,[1,2,5,6,9,10,11,13,15]),size(fys,1),[]);

fys.data = cat(2,reshape(fys(:,1:2:end,[2,6,10,11,13,15]),size(fys,1),[]),ffet(:,[1,2,9,10,11,15,16]));

figure,
subplot(211);
imagesc(ts,fs(1:2:end),log10(ys.data(:,1:2:end,1))'),axis('xy');caxis([-5,0]);colormap('jet');
subplot(212);
plot_stc(stc,1);
linkaxes(findobj(gcf,'Type','Axes'),'x');

%fys.data = reshape(cat(2,fys(:,1:3:end-16,:),sum(fys(:,1:3:end-16,:),2)),size(fys,1),[]);


%sper = stc{'a'};
sper = any(stcm.data,2);
sper([1:10,end-10:end]) = false;
[mappedX] = tsne(fys(sper,:),sum(stcm(sper,:),2),2,10,50);




figure();
hold('on');
sclr = 'rbgcmy';
for s = 1:6,
    try
    plot(mappedX(stcm(sper,s)==s,1),mappedX(stcm(sper,s)==s,2),[sclr(s),'.']);
    end
end



sper = [stc{'w+n+p'}];
[mappedWNP] = tsne(fys(sper,:),sum(stcm(sper,:),2)+1,2,10,50);


figure,
plot(mappedWNP(:,1),mappedWNP(:,2),'.');



figure();
hold('on');
sclr = 'rbgcmy';
for s = 1:6,
    try
    plot(mappedWNP(ctcm(sper,s)==s,1),mappedWNP(ctcm(sper,s)==s,2),[sclr(s),'.']);
    end
end



% plot the tsne map for all states -> examine structure.

sper = stc{'a'};
fys = ys.copy();
fys.data = reshape(cat(2,log10(fys(:,1:3:end-16,:)),sum(log10(fys(:,1:3:end-16,:)),2)),size(fys,1),[]);
[mappedXP] = tsne(fys(sper,:),sum(stcm(sper,:),2)+1,2,10,50);



sper = [stc{'w+n'}];
[mappedWN] = tsne(fys(sper,:),[],2,10,50);

% plot the tsne map for all states -> examine structure.



skeys = 'wnprms';
sclr = 'bgcrmk';

figure();
hold('on');
for s = 1:numel(skeys),
    subplot(1,numel(skeys),s);
    wys = log10(ys(stc{skeys(s)},:));
    plot(fs,wys(1:5:end,:)',sclr(s));
    ylim([-10,2]);
    title(stc{skeys(s)}.label);
end
