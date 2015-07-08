function unit_lfp_spkcohere(Trial,varargin)

%Trial = MTATrial('jg05-20120317');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;




hunits = unique(Trial.spk.create(Trial,rhm.sampleRate,'hswalk').clu)';Trial.spk.clear;
lunits = unique(Trial.spk.create(Trial,rhm.sampleRate,'lswalk').clu)';Trial.spk.clear;
runits = unique(Trial.spk.create(Trial,rhm.sampleRate,'rear').clu)';Trial.spk.clear;
Trial.spk.create(Trial,rhm.sampleRate);

%%%%
% $$$ [yt, ft,  phit, yerrt, phierrt, phloct,powt] =  mtptchd(wrhm,Trial.spk(109),ones(numel(Trial.spk(109)),1) , 2^9,rhm.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30],[],Trial.stc{'l',rhm.sampleRate}.data);%,pval,Mi
%%%%

lfp = Trial.lfp.copy;
lfp.load(Trial,65:96);
wlfp = WhitenSignal([lfp.data],[],1);

spk = {};
states = {'theta','rear','walk'};
nsts = numel(states);

for i = 1:nsts,
    spk{i} = Trial.spk.copy;
    spk{i}.create(Trial,lfp.sampleRate);

    [yl, fl,  phil, yerrl, phierrl, phlocl,powl] =  mtptchd(wrhm,spk{i}.res, Trial.spk.clu, 2^9,rhm.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30],[],Trial.stc{'l',rhm.sampleRate}.data,[],10);
    yls = zeros([size(yl,1),size(Trial.spk.map,1)]);
    ylt = zeros([size(yl,1),lfp.size(2),size(Trial.spk.map,1)]);
    philt = zeros([size(yl,1),lfp.size(2),size(Trial.spk.map,1)]);
    for i=1:numel(units)
        yls(:,units(i)) = yl(:,i+lfp.size(2),lfp.size(2));
        for j=1:lfp.size(2)
            ylt(:,j,units(i)) = yl(:,i+lfp.size(2),j);
            philt(:,j,units(i)) = phil(:,i+lfp.size(2),j);
        end
    end

    pfs{i} = MTAApfs(Trial,[],states{i},true);
end


[accg,tbin] = autoccg(Trial);

figH = figure(32391);
set(figH,'paperposition',get(figH,'position').*[0,0,1,1]./30)
for i = 1:size(Trial.spk.map,1),
    clf,
    subplot(4,3,1);cla
    pfr.plot(i);
    title(['rear unit:' num2str(i)]);
    subplot(4,3,2);cla
    pfg.plot(i);
    title(['hwalk unit:' num2str(i)]);
    subplot(4,3,3);cla
    pfl.plot(i);
    title(['lwalk unit:' num2str(i)]);
    subplot(4,3,4);cla
    imagescnan({fr,1:32,yrt(:,:,i)'},[],0,1);
    title('LFP X Unit Coherence')
    ylabel('channels')
    xlabel('frequency')
    subplot(4,3,5),cla
    imagescnan({fg,1:32,ygt(:,:,i)'},[],0,1);
    title('LFP X Unit Coherence')
    xlabel('frequency')
    subplot(4,3,6);cla
    imagescnan({fl,1:32,ylt(:,:,i)'},[],0,1);
    title('LFP X Unit Coherence')
    xlabel('frequency')
    subplot(4,3,7);cla
    imagescnan({fr,1:32,phirt(:,:,i)'},[],1,1);
    title('LFP X Unit Phase Diff')
    ylabel('channels')
    subplot(4,3,8);cla
    imagescnan({fg,1:32,phigt(:,:,i)'},[],1,1);
    title('LFP X Unit Phase Diff')
    subplot(4,3,9);cla
    imagescnan({fl,1:32,philt(:,:,i)'},[],1,1);
    title('LFP X Unit Phase Diff')
    subplot(4,3,10);cla
    bar(tbin,accg(:,i));axis tight
    subplot(4,3,11);cla
    title('RHM X Unit Coherence')
    hold on
    plot(fl,yrs(:,i),'r');
    plot(fg,ygs(:,i),'b');
    plot(fl,yls(:,i),'g');
    %legend('rear','lwalk','hwalk')
    subplot(4,3,12);
    title('LFP X Unit Coherence Pyr')
    hold on
    plot(fl,yrt(:,6,i),'r');
    plot(fg,ygt(:,6,i),'b');
    plot(fl,ylt(:,6,i),'g');
    xlabel('frequency')
    h = legend('r','h','l');
    %waitforbuttonpress
    %pause(1)

    saveas(figH,fullfile(Trial.path.data,'figures',[Trial.filebase,'-unit_rhm_spkcohere'],[Trial.filebase,'-unit_rhm_spkcohere-',num2str(i),'.png']));
    %reportfig(fullfile(Trial.path.data,'figures'),figH, ...
    %[Trial.filebase '-unit_rhm_spkcohere'],[],[Trial.filebase ' unit:' num2str(i)],300)
end



