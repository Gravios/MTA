function [sper] = bhv_shake(Trial)

% Trial = MTATrial('jg05-20120317');
% Trial = MTATrial('Ed05-20140529','all','ont');
% Trial = MTATrial('Ed01-20140709');
% Trial = MTATrial('Ed03-20140625');


xyz = Trial.load('xyz');
ix = nniz(xyz(:,1,1));
xyz.data = ButFilter(xyz.data,3,[55]/(xyz.sampleRate/2),'low');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
hcom = xyz.com(rb);
bcom = xyz.com(rbb);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('fbcom',[.7,1,.7],{{'spine_lower','spine_middle',[0,0,1]}},ButFilter(bcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
ang = create(MTADang,Trial,xyz);

sparm = struct('nFFT',2^7,...
               'Fs',ang.sampleRate,...
               'WinLength',2^5,...
               'nOverlap',2^5*.875,...
               'FreqRange',[5,50]);

fet = xyz.copy;
fet.ext = 'fet';
fet.label = 'shake';
fet.key = 'k';
fet.data = [diff(circ_dist(circ_dist(ang(:,3,11,1),ang(:,3,4,1)),pi).*double(ix)),...
            diff(circ_dist(circ_dist(ang(:,3,5,1),ang(:,3,10,1)),0).*double(ix)),...
            diff(circ_dist(circ_dist(ang(:,1,11,1),ang(:,2,4,1)),0).*double(ix))];;
[ys,fs,ts] = fet_spec(Trial,fet,...
                      'mode','mtchglong',...
                      'wsig',true,...
                      'defspec',sparm,...
                      'overwrite',false);

% load normalization vars
load('MTAFet_model_shake_k.mat'
ysMean = [];
ysStd = [];

find = fs<25&fs>10;
[ys.data,ysMmean,ysStd] = nunity(cat(2,median(log10(ys(:,find,1,1)),2),...
                                       median(log10(ys(:,find,2,2)),2),...
                                       median(log10(ys(:,find,3,3)),2)),...
                                 [],ysMean,ysStd);
        
[sfet,sMean,sStd] = nunity(prod(ys.data,2));


% fet.data = diff(circ_dist(circ_dist(ang(:,3,11,1),ang(:,3,4,1)),pi).*double(ix));
% [ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',sparm);
% 
% fet.data = diff(circ_dist(circ_dist(ang(:,3,5,1),ang(:,3,10,1)),0).*double(ix));
% [rs,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',sparm);
% 
% fet.data = diff(circ_dist(circ_dist(ang(:,1,11,1),ang(:,2,4,1)),0).*double(ix));
% [gs,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',sparm);
% 
% figure,imagesc(ts,fs,log10(gs.data)'),axis xy
% Lines(Trial.stc{'k',1}(:),[],'m');



hfig = figure(11);clf,hold on,
set(hfig,'PaperType','A3')
set(hfig,'paperposition',[0,0,10,5])
plot(ts,nunity(log10(median(ys(:,fs<20&fs>10),2)))),
plot(ts,nunity(log10(median(rs(:,fs<20&fs>10),2))),'r'),
plot(ts,nunity(log10(median(gs(:,fs<20&fs>10),2))),'c'),
Lines(Trial.stc{'k',1}(:),[],'m');
Lines(Trial.stc{'r',1}(:),[],'r');
Lines([],-5.5,'b');
Lines([],-6.5,'r');
title([Trial.filebase ': Shake Feature k_1 and k_2'])
xlabel('Time (s)')
legend('k_1','k_2')
saveas(hfig,fullfile('/storage/gravio/figures/BHV_detection/Shake',...
                     [Trial.filebase '-Fet-k_1_2.eps']),'epsc');

ind = Trial.stc{'a'};
[~,Am,As] = nunity(log10([median(ys(ind,fs<20&fs>10),2),...
                    median(rs(ind,fs<20&fs>10),2)]));

ind = Trial.stc{'a'};
figure
hist2(nunity(log10([median(ys(ind,fs<20&fs>10),2),...
                    median(rs(ind,fs<20&fs>10),2)]),[],Am,As),...
      linspace(-3,6,100),...
      linspace(-3,6,100))

ind = Trial.stc{'a'};
figure,plot(nunity(log10([median(ys(ind,fs<20&fs>10),2),...
                    median(rs(ind,fs<20&fs>10),2)]),[],Am,As))

hist2(nunity(log10([median(ys(ind,fs<20&fs>10),2),...
                    median(rs(ind,fs<20&fs>10),2)]),[],Am,As),...
      linspace(-3,6,100),...
      linspace(-3,6,100))

                

