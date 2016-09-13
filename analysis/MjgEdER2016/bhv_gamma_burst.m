


Trial = MTATrial.validate('jg05-20120310.cof.all');

ds = load(fullfile(Trial.spath,[Trial.name,'.DetectGammaBursts3.lfpinterp.all.1-96.mat']));




[~,tind] = SelectPeriods(ds.BurstTime,Trial.sync.data,'d',1,0);

for field = fieldnames(ds)'
    field = field{1};
    if ~isstruct(ds.(field))
        ds.(field) = ds.(field)(tind);
    end
end

channels = 65:96;

ds.BurstTime = ds.BurstTime-Trial.sync(1);

% $$$ s = 'p';
% $$$ channels = 65:96;
% $$$ fbins = ds.Params.SpecF;
% $$$ sper = Trial.stc{s,1};
% $$$ [~,bind] = SelectPeriods(ds.BurstTime,sper.data,'d',1,0);
% $$$ binds = false([numel(ds.BurstTime),1]);
% $$$ binds(bind)=true;
% $$$ bchanc=[];
% $$$ for i = channels,
% $$$     bchanc(:,end+1) = histc(ds.BurstFreq(ismember(ds.BurstChan,i)&binds),fbins);
% $$$ end
% $$$ figure,imagesc(fbins,channels,bchanc')



pchans = [69,73,81,96];
% Load LFP
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',pchans);
tbp_phase = lfp.phase;

flfp = lfp.copy;
flfp.filter('ButFilter',3,[5,13],'bandpass');

[minds,mvals] = LocalMinima(-flfp(:,1),50,-1000);

minds = (minds-1)/flfp.sampleRate;

[xx,xi] = NearestNeighbour(minds,ds.BurstTime);

xi(xi==1)=2;

ds.BurstThetaFreq = 1./mean(diff(GetSegs(minds,xi-1,3)))';
ds.BurstThetaPhase = tbp_phase(round(ds.BurstTime.*lfp.sampleRate)+1,:);
ds.BurstThetaEnvelope = 1./mean(diff(GetSegs(mvals,xi-1,3)))';


s = 'w';
channels = 65:96;
fbins = ds.Params.SpecF;
sper = Trial.stc{s,1};
[~,bind] = SelectPeriods(ds.BurstTime,sper.data,'d',1,0);
binds = false([numel(ds.BurstTime),1]);
binds(bind)=true;
bchanc=;
for i = channels,
    [~,bchanc(:,end+1)] = histc(ds.BurstFreq(ismember(ds.BurstChan,i)...
                                             &binds...
                                             &),fbins);

end


figure,imagesc(fbins,channels,bchanc')




xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
vxy = xyz.vel([1,5],[1,2]);
vxy.data(vxy.data<1e-3) = 1e-4;
vxy.data = log10(vxy.data);

ds.BurstHeadSpeed = vxy(round(minds(xi).*vxy.sampleRate)+1,:);

vbins = linspace(-3,2,20);
%vbinds = histc(




freqRange = [20,40;...
             30,50;...
             40,60;...
             50,70;...
             60,80;...
             70,90;80,100;90,110;100,120;110,130;120,140;130,150;140,160];
tfr = bsxfun(@plus,[1:6;2:7],[5;5])';
tph = bsxfun(@plus,[-pi:.5:pi-.25],[0;0.5])';

velRange = bsxfun(@plus,[-3:.5:1.5]',[0,0.5]);
bmp=[];
bsp=[];
bte=[];
for i = channels,
    for f = 1:size(freqRange,1),
        for r = 1:numel(pchans),
            for v = 1:size(velRange),
            %for t = 1:size(tfr,1),
                %for p = 1:size(tph,1),
                    bselind = ismember(ds.BurstChan,i)...
                              &ds.BurstFreq<freqRange(f,2)...
                              &ds.BurstFreq>freqRange(f,1)...
                              &ds.BurstHeadSpeed(:,2)<velRange(v,2)...
                              &ds.BurstHeadSpeed(:,2)>velRange(v,1);
                    %&ds.BurstThetaFreq<tfr(t,2)...
                    %&ds.BurstThetaFreq>tfr(t,1);...

                    %&ds.BurstThetaFreq<tfr(t,2)...
                    %&ds.BurstThetaFreq>tfr(t,1);...
                    %&ds.BurstThetaPhase(:,1)<tph(p,2)...
                    %          &ds.BurstThetaPhase(:,1)>tph(p,1);
                    bmp(f,i-64,r,v) = circ_mean(ds.BurstThetaPhase(bselind,r));
                    bsp(f,i-64,r,v) = circ_std(ds.BurstThetaPhase(bselind,r));
%bte(f,i-64,r,t,p) = mean(ds.BurstThetaEnvelope(bselind,1));
%bte(f,i-64,t,p) = mean(ds.BurstThetaEnvelope(bselind,1));
%end
%end
            end
        end    
    end
end


figure,
for r = 1:numel(pchans),
    for t = 1:size(velRange,1),
        subplot2(numel(pchans),size(velRange,1),r,t);
        imagesc(mean(freqRange,2),channels,bmp(:,:,r,t)');
        colormap hsv
    end
end

figure,
for r = 1:numel(pchans),
    for t = 1:size(tfr,1),
        subplot2(numel(pchans),size(tfr,1),r,t);
        imagesc(mean(freqRange,2),channels,bsp(:,:,r,t)');
        colormap jet
    end
end
ForAllSubplots('caxis([0,1.4])');


% $$$ figure,
% $$$ for r = 1:numel(pchans),
% $$$     for t = 1:size(tfr,1),
% $$$         subplot2(numel(pchans),size(tfr,1),r,t);
% $$$         imagesc(mean(freqRange,2),channels,bmp(:,:,r,t)');
% $$$         colormap hsv
% $$$     end
% $$$ end
% $$$ 
% $$$ figure,
% $$$ for r = 1:numel(pchans),
% $$$     for t = 1:size(tfr,1),
% $$$         subplot2(numel(pchans),size(tfr,1),r,t);
% $$$         imagesc(mean(freqRange,2),channels,bsp(:,:,r,t)');
% $$$         colormap jet
% $$$     end
% $$$ end
% $$$ ForAllSubplots('caxis([0,1.4])');


figure,
for r = 1;%:numel(pchans),
    for t = 1:size(tfr,1),
        subplot(1,size(tfr,1),t);
        imagesc(mean(freqRange,2),channels,bsp(:,:,r,t)');
        colormap jet
    end
end


figure,
for r = 1:numel(pchans),
subplot(1,numel(pchans),r);
    imagesc(mean(freqRange,2),channels,bsp(:,:,r)');
colormap jet
end



figure,
c = 16;
for f = 1:size(freqRange,1),
    subplot(1,size(freqRange,1),f);
    imagesc(mean(tph,2),mean(tfr,2),sq(bte(f,c,:,:))');
    colormap jet
    caxis([-.1,.1])
end








freqRange = [20,40;...
             30,50;...
             40,60;...
             50,70;...
             60,80;...
             70,90;80,100;90,110;100,120;110,130;120,140;130,150;140,160];
tfr = bsxfun(@plus,[1:6;2:7],[5;5])';
tph = bsxfun(@plus,[-pi:.5:pi-.25],[0;0.5])';

velRange = bsxfun(@plus,[-3:.5:1.5]',[0,0.5]);
bmp=[];
bsp=[];
bte=[];
for i = channels,
    for f = 1:size(freqRange,1),
        for r = 1:numel(pchans),
            for v = 1:size(velRange),
            %for t = 1:size(tfr,1),
                %for p = 1:size(tph,1),
                    bselind = ismember(ds.BurstChan,i)...
                              &ds.BurstFreq<freqRange(f,2)...
                              &ds.BurstFreq>freqRange(f,1)...
                              &ds.BurstHeadSpeed(:,2)<velRange(v,2)...
                              &ds.BurstHeadSpeed(:,2)>velRange(v,1);
                    %&ds.BurstThetaFreq<tfr(t,2)...
                    %&ds.BurstThetaFreq>tfr(t,1);...

                    %&ds.BurstThetaFreq<tfr(t,2)...
                    %&ds.BurstThetaFreq>tfr(t,1);...
                    %&ds.BurstThetaPhase(:,1)<tph(p,2)...
                    %          &ds.BurstThetaPhase(:,1)>tph(p,1);
                    bmp(f,i-64,r,v) = circ_mean(ds.BurstThetaPhase(bselind,r));
                    bsp(f,i-64,r,v) = circ_std(ds.BurstThetaPhase(bselind,r));
%bte(f,i-64,r,t,p) = mean(ds.BurstThetaEnvelope(bselind,1));
%bte(f,i-64,t,p) = mean(ds.BurstThetaEnvelope(bselind,1));
%end
%end
            end
        end    
    end
end


rper = Trial.stc{'r'};
bc = ds.BurstChan==73;
bt = round(ds.BurstTime(bc).*rper.sampleRate)+1;
bf = round(ds.BurstFreq(bc),-1);
[trl,tri,trc] = TrigRasters(rper(:,1)-60,120,bt,bf,rper.sampleRate,1);
