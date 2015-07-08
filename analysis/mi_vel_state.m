
Trial = MTATrial('jg05-20120317')
state = 'w';
Trial.load('xyz')
Trial.load('ang')

dwin= gausswin(61)./sum(gausswin(61));
fwin = gausswin(11)./sum(gausswin(11));
Trial.xyz.filter(fwin);

fpv = cat(1,[-2,-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','spine_upper','head_front'},[1,2]))));
fpz = cat(1,[-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','head_back'},[3]))));

fet =[];
fet = [fet,fpv,fpz];
fet = [fet,Trial.xyz(:,7,3)-Trial.xyz(:,1,3)];
fet = [fet,Trial.ang(:,3,4,2)];
fet = [fet,Trial.ang(:,5,7,2)];
fet = [fet,fet.^2,fet.^3];
% $$$ fet = MTADxyz([],[],Filter0(fwin,fet),Trial.xyz.sampleRate);
% $$$ fet.data(fet.data==0) = nan;
% $$$ fet_not_nan = ~isnan(fet.data(:,1));



for i = 1:size(fet,2)
    %for i = 1:fet.size(2)
    edges(:,i) = linspace(prctile(fet(~isinf(fet(:,i)),i),[.1]),prctile(fet(~isinf(fet(:,i)),i),[99.9]),32);
end

sbound = -60:60;
%ixy = repmat({zeros(numel(sbound),fet.size(2))},[fet.size(2),1]);
ixy = repmat({zeros(numel(sbound),size(fet,2))},[size(fet,2),1]);

if matlabpool('size')==0,matlabpool open 12,end

% $$$ parfor m = 1:fet.size(2)
% $$$     for o = 1:fet.size(2)
parfor m = 1:size(fet,2)
    for o = 1:size(fet,2)
        for s = 1:numel(sbound)

            vx = MTADxyz('data',fet(:,m),'sampleRate',Trial.xyz.sampleRate);
            sx = vx(Trial.stc{state});
            vy = MTADxyz('data',circshift(fet(:,o),sbound(s)),'sampleRate',Trial.xyz.sampleRate);
            sy = vy(Trial.stc{state});

            noinf = ~isinf(sy)&~isinf(sx);
            nind = sum(noinf);

            [out,xb,yb,p]=hist2([sx(noinf),sy(noinf)],edges(:,m),edges(:,o));
            pxy = out./nind;

            px = histc(sx(~isinf(sx)&~isinf(sy)),xb);
            px = px(1:end-1)/nind;

            py = histc(sx(noinf),yb);
            py = py(1:end-1)/nind;

            pixy = pxy.*log2(pxy./(px*py'));
            pixy(isinf(pixy))=nan;
            ixy{m}(s,o) = nansum(pixy(:));

        end
    end
end




nixy = reshape(cell2mat(ixy),[numel(sbound),size(fet,2),size(fet,2)]);

[mixy,sixy] = max(nixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);


figure,
subplot2(1,3,1,1),imagesc(sq(nixy(61,:,:)))
subplot2(1,3,1,2),imagesc(mixy(:,:))
subplot2(1,3,1,3),imagesc(sixy(:,:)./Trial.xyz.sampleRate.*1000)



figure
subplot(3,1,1),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,7,:))'),axis xy
subplot(3,1,2),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,4,:))'),axis xy
subplot(3,1,3),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,1,:))'),axis xy

%%log10 version
figure
subplot(3,1,1),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,7,:))')),axis xy
subplot(3,1,2),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,4,:))')),axis xy
subplot(3,1,3),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,1,:))')),axis xy



rp = Trial.stc{state}.copy;rp.cast('TimeSeries');
wp = Trial.stc{'w'}.copy;wp.cast('TimeSeries');

f1 = 24;
f2 = 22;
nif = ~isinf(fet(:,f1))&~isinf(fet(:,f2))&rp.data;
nif = ~isinf(fet(:,f1))&~isinf(fet(:,f2))&wp.data;
nif = ~isinf(fet(:,f1))&~isinf(fet(:,f2));
figure,hist2([clip(fet(nif,f1),edges(1,f1),edges(end,f1)),clip(fet(nif,f2),edges(1,f2),edges(end,f2))],50,50)

%figure,hist2([10.^clip(fet(nif,f1),edges(1,f1),edges(end,f1)),10.^clip(fet(nif,f2),edges(1,f2),edges(end,f2))],50,50)


fet(isinf(fet))=nan;