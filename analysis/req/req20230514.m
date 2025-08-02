

MjgER2016_load_data();



Trial = Trials{20};

xyz = preproc_xyz(Trial,'trb');
vxy = vel(filter(copy(xyz),'ButFilter',4,2,'low'),'hcom',[1,2]);
lvxy = vel(filter(copy(xyz),'ButFilter',4,0.1,'low'),'hcom',[1,2]);


stc = Trial.stc.copy();
units = select_units(Trial,'int');
spk = Trial.load('spk',xyz.sampleRate,'',units);

pper = [stc{'s',xyz.sampleRate}]

vper = ThreshCross(vxy.data,1,1);
vper(lvxy(round(mean(vper,2)))>0.4,:) = [];
%vper(vper(1:end-1,1)-vper(2:end,2)-
vper = vper(WithinRanges(round(mean(vper,2)),pper.data),:);
whos vper

%vper(diff(vper,1,2)>100,:) = [];


sper = stc{'x'};
% $$$ rper = stc{'r'};
% $$$ 

% $$$ 
% $$$ 
% $$$ per = stc.get_state_transitions(Trial,{'x','p'},0.2,xyz);
% $$$ dper = stc.get_state_transitions(Trial,{'p','x'},0.2,xyz);
% $$$ 
% $$$ 
% $$$ 
% $$$ mper = sper.data;
% $$$ mper(diff(mper,1,2)<60,:) = [];
% $$$ 
% $$$ mper = rper.data;
% $$$ mper(diff(mper,1,2)<60,:) = [];
% $$$ 
% $$$ mper = per;
% $$$ mper = dper;

mper = vper(1:2:end,:);
%mper = sper(1:2:end,:);
disp(size(vper))

figure
for u = 1:numel(units),
mspk = spk(units(u));
[mccg, bins] = CCG([mspk;mper(:,1);mper(:,2)],...
                   [ones(size(mspk)); 2*ones([length(mper),1]); 3*ones([length(mper),1])],...
                   4,...
                   40,...
                   spk.sampleRate,...
                   [1,2,3],...
                   'scale');
subplot(5,8,u);
bar(bins,mccg(:,1,3));
title(['Mov On ',num2str(units(u))])
ylim([0,2.5]);
end

figure
for u = 1:numel(units),
mspk = spk(units(u));
[mccg, bins] = CCG([mspk;mper(:,1);mper(:,2)],...
                   [ones(size(mspk)); 2*ones([length(mper),1]); 3*ones([length(mper),1])],...
                   4,...
                   40,...
                   spk.sampleRate,...
                   [1,2,3],...
                   'scale');
subplot(5,8,u);
bar(bins,mccg(:,1,1));
title(['Mov Off ',num2str(units(u))])
ylim([0,2.5]);
end

figure
subplot(211);
hold('on');
plot([1:size(vxy,1)]./vxy.sampleRate,(vxy.data))
plot([1:size(vxy,1)]./vxy.sampleRate,(lvxy.data))
Lines(mper(:)./vxy.sampleRate,[],'k');
subplot(212);
plotSTC(stc,1);
linkx();


rper = stc{'R-t&s',vxy.sampleRate};
mper = rper.data(1:3:end,:);

figure
for u = 1:numel(units),
mspk = spk(units(u));
[mccg, bins] = CCG([mspk;mean(mper(:,[1,2]),2)],...
                   [ones(size(mspk)); 2*ones([length(mper),1]); ],...
                   1,...
                   80,...
                   spk.sampleRate,...
                   [1,2],...
                   'scale');
subplot(5,8,u);
bar(bins,mccg(:,1,2));
title(['SRW ',num2str(units(u))])
ylim([0,2.5]);
end


remPer = [stc{'s&t',xyz.sampleRate}];   label = 'REM';
remPer = [stc{'t-s-m',xyz.sampleRate}]; label = 'RUN';
remPer = [stc{'s-t',xyz.sampleRate}];   label = 'REST';
remPer.data(:,2) = remPer.data(:,1)+400;
remPer.data(:,1) = remPer.data(:,1)+100;

ind = 1:3;
%ind = ':';
figure
for u = 1:numel(units),
mspk = spk(units(u));
mspk = mspk(WithinRanges(mspk,remPer.data(ind,:)));
[mccg, bins] = CCG([mspk;remPer(ind,1);remPer(ind,2)],...
                   [ones(size(mspk)); 2*ones([size(remPer(ind,1),1),1]); 3*ones([size(remPer(ind,2),1),1])],...
                   1,...
                   80,...
                   spk.sampleRate,...
                   [1,2,3],...
                   'hz');
subplot(5,8,u);
bar(bins,mccg(:,1,1));
title([label,' ',num2str(units(u))])
%ylim([0,2.5]);
end



u = 5;
accg = zeros([size(xyz,1),193,numel(units)]);
for u = 1:numel(units)
mspk = spk(units(u));

trainNextIndex = 1;
for ind = 129:32:size(xyz,1)-129,
    twin = [-128,128]+ind;
% $$$     trainStartIndex = trainNextIndex;
% $$$     trainEndIndex = trainStartIndex;
% $$$     while true
% $$$         if mspk(trainEndIndex) < twin(2)
% $$$             trainEndIndex = trainEndIndex + 1;
% $$$         else
% $$$             trainNextIndex = trainEndIndex;
% $$$             break
% $$$         end
% $$$     end
% $$$     if (trainEndIndex - trainStartIndex) < 2 
% $$$         continue
% $$$     end
% $$$     tspk = mspk(trainStartIndex:trainEndIndex);
    tspk = mspk(WithinRanges(mspk,twin));
    [mccg, bins] = CCG([tspk],...
                       [ones(size(tspk))],...
                       1,...
                       96,...
                       spk.sampleRate,...
                       [1],...
                   'hz');
    accg(ind,:,u) = mccg(:,1,1)';
end
end

figure,imagesc(accg(1:32:end,:)')
    

% $$$ FreqRange = [0,xyz.sampleRate/2];
% $$$ fo = ([1:256/2+1]-1)'*xyz.sampleRate/256;
% $$$ select = find( fo > FreqRange(1,1) & fo < FreqRange(1,2));
% $$$ fo = fo(select);


fftOutAll = nan([size(xyz,1),127,numel(units)]);
for u = 1:numel(units)
for i = 129:32:size(xyz,1)-129
    fftOut = fft(accg(i,:,u)',256).*sqrt(2);
    fftOutAll(i,:,u) = fftOut(select).*conj(fftOut(select));
end
end

figure,plot(imgaussfilt(log10(fftOutAll(1:32:end,1,u))))

u = 6
figure,imagesc((1:32:size(fftOutAll,1))./xyz.sampleRate,fo,...
               bsxfun(@rdivide,...
                      imgaussfilt(log10(fftOutAll(1:32:end,:,u))',1.5),...
                      sum(imgaussfilt(log10(fftOutAll(1:32:end,:,u))',1.5))...
                      ) ...
               );

figure,imagesc((1:32:size(fftOutAll,1))./xyz.sampleRate,fo,...
               bsxfun(@rdivide,...
                      imgaussfilt(log10(mean(fftOutAll(1:32:end,:,[5,6,11,13]),3,'omitnan'))',1.5),...
                      sum(imgaussfilt(log10(mean(fftOutAll(1:32:end,:,[5,6,11,13]),3,'omitnan'))',1.5))...
                      ) ...
               );
axis('xy');
colormap('jet');
title(num2str(units(u)))
% $$$ figure,imagesc((1:64:size(fftOutAll,1))./xyz.sampleRate,fo,bsxfun(@rdivide,imgaussfilt(log10(fftOutAll(1:64:end,:))',2),sum(imgaussfilt(log10(fftOutAll(1:64:end,:))',2))));
%caxis([3,8]);   
caxis([0.007,0.011])
colormap('jet');
tper = [stc{'t',1}]
Lines(tper.data(:,1),[],'g');
Lines(tper.data(:,2),[],'m');

accg15 = accg;

figure,imagesc((1:60:size(fftOutAll,1))./xyz.sampleRate,...
               fo,...
               bsxfun(@rdivide,log10(fftOutAll(1:60:end,:))'
axis('xy');
%caxis([3,8]);   
caxis([0.006,0.013])
colormap('jet');




size(CCG([mspk],...
                       [ones(size(mspk))],...
                       1,...
                       80,...
                       spk.sampleRate,...
                       [1],...
                   'hz'))



u = 13;
figure,
subplot (221) ,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,18,u))],linspace([-3,2,50]),linspace([-1,8,50]),'xprob'),caxis([0,0.1,]),colormap('jet');
subplot (222) ,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,1,u))],linspace([-3,2,50]),linspace([2,8,50]),'yprob'),caxis([0,0.1,]),colormap('jet');
subplot (223) ,hist2([log10(fftOutAll(1:32:end,1,u)),log10(fftOutAll(1:32:end,14,u))],linspace([3,8,50]),linspace([-1,8,50]),'yprob'),caxis([0,0.1,]),colormap('jet');



figure,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,18,u))],linspace([-3,2,50]),linspace([-1,8,50]))



figure,plot(log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,1,u)),'.');

figure,plot(log10(vxy(1:32:end)),sum(log10(fftOutAll(1:32:end,30:58,u)),2),'.');


figure,plot(log10(fftOutAll(1:32:end,12,u)),log10(fftOutAll(1:32:end,18,u)),'.');
