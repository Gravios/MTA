
lfp = load('/storage/weiwei/matlab/EMG_removing/EMG_rm/data_drafts/lfp1.mat');

lfp = load('/storage/weiwei/matlab/EMG_removing/EMG_rm/data_drafts/lfp1.mat');

%  Example of EMG noise detection and removing: 
%% Data Loading
load lfp1.mat %  this is the first continuse period of LFP recording, at


[x, Ws, As, EMG_au] = EMG_rm(lfp.data, lfp.sr);



%% EMG Priod Selection
tmp_t = 9e6:9.5e6; % example period with high muscle tone
[x, Ws, As, EMG_au] = EMG_rm(lfp.data, lfp.sr);
figure;clf
opf = @(x)(bsxfun(@plus, x/2e4, 1:size(x,2)));
plot(tmp_t, opf(lfp.data(tmp_t,:)), 'k'); 
hold on; 
plot(tmp_t, opf(x(tmp_t,:)), 'r');
 

EMG_rm_main('jg05-20120312',[],9);

return
%% Validation with Ach, given the Ach data.
% Ach Signals and movement related stuffs
% plot of spectrum, Ach recording and ripple power: 
load ACh_NREM.mat
opf1 = @(x)(bsxfun(@rdivide,x, std(x))); % not perticularily important.
opf = @(x)(opf1(bsxfun(@minus, x, min(x))));
end_t = length(lfp.data)/1000*5;
figure;
plot([1:end_t]/5, bsxfun(@plus, [10 20], opf([ACh_NREM.signals.ACh(1:end_t), ACh_NREM.signals.SWRpower(1:end_t)])));
hold on
plot([1:10:length(lfp.data)]/1000, opf(EMG_au(1:10:end)));
plot(ACh_NREM.SWRs/ACh_NREM.lfpSampRate, zeros(size(ACh_NREM.SWRs))-20,'r+')
plot(ACh_NREM.peaksACh/Achfs, zeros(size(ACh_NREM.peaksACh))-25,'k+')
legend('ripple peak power', 'Ach','IC', 'Ach peak', 'ripple peak')
xlim([0 1.4e4])

%filename = '/storage/gravio/data/processed/nlx/jg05-20120312/jg05-20120312.dat';
filename = '/storage/gravio/data/processed/nlx/jg05-20120312/jg05-20120312.lfp';
Par = LoadPar('/storage/gravio/data/processed/nlx/jg05-20120312/jg05-20120312.xml');
lfp = LoadBinary( filename,65:Par.nChannels,Par.nChannels,4)';

filename = '/storage/gravio/data/project/general/jg05-20120312/jg05-20120312.lfp';
x =  LoadBinary( [filename,'d'],65:Par.nChannels,Par.nChannels,4)';
% $$$ addpath(genpath('/storage/weiwei/matlab/EMG_removing/'));
% $$$ [x, Ws, As, EMG_au] = EMG_rm(lfp(:,65:96), Par.lfpSampleRate);


figure();
    hold('on');
    PlotTraces(lfp(1:1e6,1:2:32),[1:1e6]/1250,[],5,'k');
    PlotTraces(x(1:1e6,1:2:32),[1:1e6]/1250,[],5,'r');

        

    
    


ind = ':';
yll = [];
[~,arm] = WhitenSignal([x(ind,1)],[],1,[],5);
for c = 31:32
[yll(:,:,c),fll,tll] = mtcsdglong(WhitenSignal([x(ind,c)],[],1,arm,5),...
                        2^13,...
                        1250,...
                        2^12,...
                        2^12.*0.875,...
                        [],...
                        [],...
                        2,...
                        [0.1,15]);
end

figure,
imagesc(tll,fll,log10(yll(:,:,15))');
caxis([0,1.8])
colormap('jet')


figure,
subplot(411);
imagesc(tll,1:32,sq(mean(log10(yll(:,5:9,:)),2))');
caxis([0,1.8])
colormap('jet')
subplot(412);
imagesc(tll,1:32,sq(mean(log10(yll(:,20:30,:)),2))');
caxis([0,1.8])
colormap('jet')
subplot(413);
imagesc(tll,1:32,sq(mean(log10(yll(:,45:65,:)),2))');
caxis([0.5,2.5])
colormap('jet')
subplot(414);
imagesc(tll,1:32,sq(mean(log10(yll(:,45:65,:)),2))'./sq(mean(log10(yll(:,[1:40,70:98],:)),2))');
colormap('jet')
linkaxes();

ind = 15e6:15.1e6;
ysl = [];
[~,arm] = WhitenSignal([x(ind,1)],[],1,[],5);
for c = 1:32
[ysl(:,:,c),fsl,tsl] = mtcsdglong(WhitenSignal([x(ind,c)],[],1,arm,5),...
                        2^8,...
                        1250,...
                        2^7,...
                        2^7.*0.875,...
                        [],...
                        [],...
                        2,...
                        [40,280]);
end

% $$$ 
% $$$ [ym,fm,tm] = mtcsdglong(WhitenSignal([lfp(ind,82),x(ind,16)],[],1),...
% $$$                         2^9,...
% $$$                         1250,...
% $$$                         2^8,...
% $$$                         2^8.*0.875,...
% $$$                         [],...
% $$$                         [],...
% $$$                         2,...
% $$$                         [15,35]);

ind = 12570000:26385098;
yl = mtcsdglong(WhitenSignal(x(ind,c)),...
                        2^11,...
                        1250,...
                        2^10,...
                        2^10.*0.875,...
                        [],...
                        [],...
                        [],...
                        [1,35]);
yl = cat(3,yl,zeros([size(yl,1),size(yl,2),31]));
for c = 2:32,
    disp(['c: ',num2str(c)]);
    tic
    [yl(:,:,c),fl,tl] = mtcsdglong(WhitenSignal(x(ind,c)),...
                                   2^11,...
                                   1250,...
                                   2^10,...
                                   2^10.*0.875,...
                                   [],...
                                   [],...
                                   [],...
                                   [1,35]);
    toc
end


figure();
imagesc(tl,1:32,sq(mean(log10(yl(:,5:10,:)),2))');
colormap('jet');


figure();

figure();
subplot(311);
imagesc(tl,1:32,sq(mean(log10(yl(:,[1:30],:)),2))');
colormap('jet');
subplot(312);
imagesc(tl,1:32,sq(mean(log10(yl(:,6:14,:)),2))');
colormap('jet');
subplot(313);
imagesc(tl,1:32,(sq(mean(log10(yl(:,6:14,:)),2))./sq(mean(log10(yl(:,[1:4,20:30],:)),2)))');
colormap('jet');
caxis([0.5,2]);
linkaxes();


figure();
subplot(211);
imagesc(tl,fl,log10(yl(:,:,16))');


sp = [];
figure,
sp(1) = subplot(311);
imagesc(ts,fs,log10(ys(:,:,1,1))');
sp(2) = subplot(312);
imagesc(ts,fs,log10(ys(:,:,2,2))');
sp(3) = subplot(313);
hold('on');
PlotTraces(lfp(1:1e6,65:2:96),[1:1e6]/1250,[],5);
PlotTraces(x(1:1e6,1:2:32),[1:1e6]/1250,[],5,'r');
linkaxes(sp,'x');
colormap('jet');


sp(3) = subplot(322);
hold('on');
PlotTraces(lfp(1:1e6,65:2:96),[1:1e6]/1250,[],5);
PlotTraces(x(1:1e6,1:2:32),[1:1e6]/1250,[],5,'r');
linkaxes(sp,'x');
colormap('jet');

ts = ts+diff(ts(1:2))/2;
tl = tl+diff(tl(1:2))/2;

figure,
sp=[];
sp(end+1)=subplot(711);
imagesc(ts,fs,log10(sq(ysl(:,1,:)))');  title(num2str(fs(1)));
sp(end+1)=subplot(712);
imagesc(ts,fs,log10(sq(ysl(:,5,:)))');  title(num2str(fs(5)));
sp(end+1)=subplot(713);
imagesc(ts,fs,log10(sq(ysl(:,10,:)))');  title(num2str(fs(10)));
sp(end+1)=subplot(714);
imagesc(ts,fs,log10(sq(ysl(:,20,:)))');  title(num2str(fs(20)));
sp(end+1)=subplot(715);
imagesc(ts,fs,log10(sq(ysl(:,30,:)))');  title(num2str(fs(30)));
sp(end+1)=subplot(716);
imagesc(tl,fl,log10(sq(yl(:,8,:))'));        title('theta')
sp(end+1)=subplot(717);
PlotTraces(x(ind-ind(1)+1,1:2:32),[ind-ind(1)+1]/1250,[],5,'r');
linkaxes(sp,'x');
colormap('jet');
ForAllSubplots('caxis([2,3.5]);')
caxis(sp(6),[1.5,2.8])
%caxis(sp(6),[-1.5,1.5])
%xlim([25,30])

figure,
imagesc(ts,fs,log10(sq(ysl(:,30,1,1,:)))');  title(num2str(fs(30)));
xlim([       25,30])
colormap('jet');
caxis([2,4])

sp(3) = subplot(322);
hold('on');
PlotTraces(lfp(1:1e6,65:2:96),[1:1e6]/1250,[],5);
PlotTraces(x(1:1e6,1:2:32),[1:1e6]/1250,[],5,'r');
linkaxes(sp,'x');
colormap('jet');

