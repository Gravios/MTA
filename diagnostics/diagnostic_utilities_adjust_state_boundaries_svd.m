

% TIME vectors
wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);



i = 1;
figure();
hold('on');
plot(fetSegs(:,:,i),'b');
plot(nanmean(fetSegs(:,:,i),2),'r');


sp = [];
figure();
sp(end+1)=subplot2(6,1,[1:4],1);
hold('on');
plot(ts{1},nansum(csw{1}(2,:,1),3)),
%plot(ts{1},nansum(csw{1}(2,:,2),3)),
%plot(ts{1},nansum(csw{1}(2,:,:),3)),
Lines(mean(sts{1},2)./sampleRate{1},[],'g');
sp(end+1)=subplot2(6,1,5,1);
plotSTC(Stc{1},1);
sp(end+1)=subplot2(6,1,6,1);
plotSTC(StcNN{1},1);
linkaxes(sp,'x');





[sccg,txx,pxx] = cf(@(s,n,sr) CCG([s;n],[ones(size(s));2*ones(size(n))],...
                                  2,40,sr,[1,2],'count'),...
                    ssmins,nsmins,sampleRate);
accg = sum(cat(4,sccg{:}),4);
medianCorrectionOffset = median(cat(1,nsmins{:})-cat(1,ssmins{:}))./sampleRate{1};

% SUPFIG Rear
figure();
bar(txx{1},accg(:,1,2));
Lines(medianCorrectionOffset.*1000,[],'r');
title(['CCG between ' subsequentState ' onsets and local minima of onset regression'])
xlabel('Time shift(ms) centered on local minima')
ylabel('count')

FigName = 'State_transition_regression_realignment_PC1_all_to_rearOn';
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));




csegs = cf(@(a,m) a.segs(m-120,360),sfet,nsmins);
figure,
corSegs = cat(2,csegs{:});
sp = [];
nHalf = round(size(fetSegs,2)/2);
fetSegsSE = zeros([size(fetSegs,1),size(fetSegs,3)]);
corSegsSE = zeros([size(corSegs,1),size(corSegs,3)]);
for i = 1:2,
sp(end+1)=subplot2(2,2,i,1);
hold('on');
plot(fetSegs(:,:,i),'b')



figure,sp=[];
sp(end+1) = subplot2(5,1,1:3,1);
plot(csw{1}(2,:,1)./5000),hold on,plot(fet{1}(:,15)),Lines(Stc{1}{'r'}(:,1),[],'r');
sp(end+1) = subplot2(5,1,4,1);
plotSTC(Stc{1});
sp(end+1) = subplot2(5,1,5,1);
plotSTC(StcCor{1});
linkaxes(sp,'xy');



sp = [];
figure();
sp(end+1)=subplot2(3,1,1,1);
plotSTC(StcHL,1);
sp(end+1)=subplot2(3,1,2,1);
plotSTC(ds.stc{s},1);
sp(end+1)=subplot2(3,1,3,1);
plotSTC(StcCor,1);
linkaxes(sp,'x');
