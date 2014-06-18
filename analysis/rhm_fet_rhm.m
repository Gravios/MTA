
%fwin = gausswin(61)./sum(gausswin(61));

%Trial = MTATrial('Ed05-20140528');
Trial = MTATrial('jg05-20120317');


fwin = gausswin(91)./sum(gausswin(91));

Trial.load('xyz');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root'});

hcom = Trial.com(rb);
Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(fwin,hcom),[1,3,2]));


bcom = Trial.com(rbb);
Trial.addMarker(Trial.xyz,'bcom',[.7,0,.7],{{'spine_lower','pelvis_root',[0,0,1]}},bcom);
Trial.addMarker(Trial.xyz,'fbcom',[.7,0,.7],{{'spine_lower','pelvis_root',[0,0,1]}},permute(Filter0(fwin,bcom),[1,3,2]));



ang = Trial.ang.copy;
ang.create(Trial);
ang = ang(:,5,11,3);
bang = ButFilter(ang,3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
dbang = ButFilter(bang,3,[5,13]./(Trial.ang.sampleRate/2),'bandpass');


Trial.lfp.load(Trial,1:2);
Trial.lfp.resample(Trial.xyz);

% $$$ Trial.lfp.load(Trial,65:96);
% $$$ phase = Trial.lfp.phase;
% $$$ phase.resample(Trial.xyz);
% $$$ 
% $$$ Trial.spk.create(Trial,Trial.xyz.sampleRate,'lwalk',[],'deburst');
% $$$ 
% $$$ rhm = Trial.lfp.copy;
% $$$ rhm.filename = '';
% $$$ rhm.sampleRate = Trial.xyz.sampleRate;
% $$$ rhmp = rhm.phase;
% $$$ 
% $$$ for i =1:106,try,hist(rhmp(Trial.spk(i)),30),end,title(num2str(i)),waitforbuttonpress,end

%wang = WhitenSignal([ang,s.lfp.data],[],1);
%wang = [ang,Trial.lfp.data];
wang = [bang,Trial.lfp.data];

%[ya,fa,ta] = mtchglong(diff(Filter0(gausswin(11)./sum(gausswin(11)),hang)),2^9,ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);
%[ya,fa,ta] = mtchglong(WhitenSignal(Filter0(gausswin(11)./sum(gausswin(11)),hang)),2^9,ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);

[ys,fs,ts,] = mtcsdglong(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);
ys = MTADlfp('data',ys,'sampleRate',1/diff(ts(1:2)));

% Speed of the head
vh = Trial.vel(11,[1,2]);
vh = MTADxyz('data',vh,'sampleRate',Trial.xyz.sampleRate);
vh.resample(ys);
edges = [0 1 2 4  8 12 16 20;...
         1 2 4 8 12 16 20 24];

edges_labels = mat2cell(edges,2,ones(1,size(edges,2)))';
chan_labels = {'Rhythmic Head Motion','Nasal Epithelium Signal','Nasal Cavity Pressure'};

vsc = [];
for i = edges,
    ind = i(1)>=vh.data&vh.data<i(2);% & 5>=mean(log10(ys(:,fs>13&fs>5,3,3)),2);
    for j = 1:size(ys,3),
        for k = 1:size(ys,3),
            if sum(ind)~=0,
                if k~=j,
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(ys(ind,:,k,j)))./mean(sqrt(ys(ind,:,k,k).*ys(ind,:,j,j)));
                else
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(ys(ind,:,k,j));
                end

            else
                vsc(find(i(1)==edges(1,:)),:) = zeros([1,size(ys,2),size(ys,3)]);
            end

        end
    end
end

figure,
for j = 1:size(ys,3),
    for k = 1:size(ys,3),
        subplot2(size(ys,3),size(ys,3),j,k),
        if j~=k
            imagesc(1:size(edges,2),fs,vsc(:,:,j,k)'),axis xy,
            title([chan_labels{j} ' & ' chan_labels{k} ' coherence'])
        else
            imagesc(1:size(edges,2),fs,log10(vsc(:,:,j,k))'),axis xy,
            title([chan_labels{j} ' PSD Binned by Head Speed'])
        end
        set(gca,'XTickLabel',cellfun(@num2str,cellfun(@transpose,edges_labels(2:2:8),'UniformOutput',false),'UniformOutput',false)');
        xlabel('Binned Head Speed cm/s')
        ylabel('Frequency Hz')
        colorbar
    end
end



yp=ys.copy;
yp.resample(Trial.xyz);
phs = Trial.lfp.phase([5,13]);
pang = LocalMinima(dbang,5);
myp = mean(log10(yp(:,fs<13&fs>5,3,3)),2);

pedges = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,;...
         0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];

apr = [];
for i = pedges,
    ind = find(i(1)<=myp & myp<i(2) & ~isnan(myp) &~isinf(myp) & myp~=0);
    if sum(ind)~=0,
        apr(find(i(1)==pedges(1,:)),:) = hist(phs(pang(ismember(pang,ind)),2),linspace(-pi,pi,15));
        apr(find(i(1)==pedges(1,:)),:) = apr(find(i(1)==pedges(1,:)),:)./sum(apr(find(i(1)==pedges(1,:)),:));
    else
        apr(find(i(1)==pedges(1,:)),:) = zeros([1,15]);
    end
end



figure,
imagesc(linspace(-pi,pi,15),pedges(2,:),apr),axis xy
caxis([0,0.12])
title('Filtered[5-13] RHM peak NCP phase Distrb')
xlabel('Phase (radians)')
ylabel('log10 NCP 5-13hz Power')





phs = Trial.lfp.phase([5,13]);
pang = LocalMinima(dbang,5);
myp = mean(log10(yp(:,fs<13&fs>5,3,3)),2);
vhp = Trial.vel(11,[1,2]);
vhp = MTADxyz('data',[0;vhp],'sampleRate',Trial.xyz.sampleRate);
 

pedges = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5;...
         0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];
edges = [0 1 2 4  8 12 16 20;...
         1 2 4 8 12 16 20 24];
vpr = [];
vps = [];
for i = pedges,
for j = edges,
    ind = find(i(1)<=myp & myp<i(2) & ~isnan(myp) & ~isinf(myp) & myp~=0 & j(1)>=vhp.data&vhp.data<j(2) );
    if sum(ind)~=0&&sum(ismember(pang,ind))~=0,
        vpr(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:) = circ_mean(phs(pang(ismember(pang,ind)),2));
        [~,vps(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:)] = circ_mtest(phs(pang(ismember(pang,ind)),2),vpr(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:));
    else
        vpr(find(i(1)==pedges(1,:)),:) = zeros([1,1,1]);
    end
end
end


figure,
imagescnan({edges(2,:),pedges(2,:),vpr},[-pi,pi],1,1,[0,0,0]),axis xy
caxis([0,0.1])
title('Filtered[5-13] RHM peak NCP phase Distrb')
xlabel('Phase (radians)')
ylabel('log10 NCP 5-13hz Power')



[ya,fa,ta,pha,fsa] = mtchglong(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],[],[],[1,30]);

figure,hold on
sp = [];
for i = 1:size(ya,3),
sp(i) = subplot(3,1,i);
imagesc(ta,fa,log10(ya(:,:,i,i)'))
axis xy
end
linkaxes(sp,'xy');


[yo,fo, phi]= mtchd(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);

figure,imagesc(ta,fa,log10(ya(:,:,1,3)')),axis xy





