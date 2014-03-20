MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317')
state = 'r';

v = Trial.vel(1:8);
v = Filter0(gausswin(61)./sum(gausswin(61)),v);

edges = linspace(-2,2,32);
sbound = -60:60;
ixy = repmat({zeros(numel(sbound),size(v,2))},[size(v,2),1]);


if matlabpool('size')==0,matlabpool open 12,end

parfor m = 1:size(v,2)
    for o = 1:size(v,2)
        for s = 1:numel(sbound)

            vx = MTADxyz('data',log10(v(:,m)),'sampleRate',Trial.xyz.sampleRate);
            sx = vx(Trial.stc{state});
            vy = MTADxyz('data',circshift(log10(v(:,o)),sbound(s)),'sampleRate',Trial.xyz.sampleRate);
            sy = vy(Trial.stc{state});

            noinf = ~isinf(sy)&~isinf(sx);
            nind = sum(noinf);

            [out,xb,yb,p]=hist2([sx(noinf),sy(noinf)],edges,edges);
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


nixy = reshape(cell2mat(ixy),[121,8,8]);
%nixy = reshape(cell2mat(ixy),[9,121,9]);
%nixy = permute(reshape(cell2mat(ixy),[9,121,9]),[2,1,3]);

[mixy,sixy] = max(nixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);


figure,
subplot2(1,3,1,1),imagesc(sq(nixy(61,:,:)))
subplot2(1,3,1,2),imagesc(mixy(:,:,f))
subplot2(1,3,1,3),imagesc(sixy(:,:,f)./Trial.xyz.sampleRate.*1000)



figure
subplot(3,1,1),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,7,:))'),axis xy
subplot(3,1,2),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,4,:))'),axis xy
subplot(3,1,3),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,sq(nixy(:,1,:))'),axis xy

%%log10 version
figure
subplot(3,1,1),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,7,:))')),axis xy
subplot(3,1,2),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,4,:))')),axis xy
subplot(3,1,3),imagesc(sbound./Trial.xyz.sampleRate.*1000,1:8,log10(sq(nixy(:,1,:))')),axis xy




