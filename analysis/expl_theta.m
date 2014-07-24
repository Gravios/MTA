Trial = MTATrial('jg05-20120310');

lfp = Trial.lfp.copy;
lfp.load(Trial,65:96);

[Data]= LoadBinary(fullfile(Trial.spath,[Trial.name '.lfp']), 65:96,96)';


wlfp = WhitenSignal(lfp.data,[],1);
wlfp = WhitenSignal(Data,[],1);


x = 0;
ys = [];
for i = 1:size(wlfp,2);
[ys(:,:,i),fs,ts] = mtchglong(wlfp((x+1):(x+1000000),i),2^10,lfp.sampleRate,2^6,2^6-2^3,3,'linear',[],[80,450]);
end

ts = ts+2^5/1250;



% $$$ figure,imagesc(ts,fs,log10(ys)'),axis xy
% $$$ nys = ys./repmat(sum(ys,2),[1,size(ys,2)]);
% $$$ figure,imagesc(ts,fs,log10(nys)'),axis xy



ndims = 2;
vect = repmat({-11:11},[ndims,1]);
sind = cell(ndims,1);
[sind{:}] = ndgrid(vect{:});
for i = 1:ndims,
    sind{i} = sind{i}.^2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));

Sys = convn(ys,Smoother,'same');
Sys = convn(Sys,Smoother,'same');
Sys = convn(Sys,Smoother,'same');

rSys = convn(Sys,Smoother,'same');
rSys = convn(rSys,Smoother,'same');
rSys = convn(rSys,Smoother,'same');
rSys = convn(rSys,Smoother,'same');
rSys = convn(rSys,Smoother,'same');
rSys = convn(rSys,Smoother,'same');


figure,sp=[];
sp(end+1)=subplot(311);imagesc(ts,fs,log10(Sys)'),axis xy
sp(end+1)=subplot(312);imagesc(ts,fs,log10(ys)'),axis xy
sp(end+1)=subplot(313);imagesc(ts,fs,log10(Sys)'-log10(rSys)'),axis xy
linkaxes(sp,'xy')
caxis([.01,.08])
xlim([181,


figure,sp=[];
sp(end+1)=subplot(311);imagesc(ts,1:32,unity(log10(sq(Sys(:,70,:))'))),%caxis([-1,1])
sp(end+1)=subplot(312);imagesc(ts,1:32,unity(log10(sq(Sys(:,80,:))'))),%caxis([-1,1])
sp(end+1)=subplot(313);imagesc(ts,1:32,unity(log10(sq(Sys(:,90,:))'))),%caxis([-1,1])
linkaxes(sp,'xy'),hold on
plot(unity(Data((1+x):(x+1000000),1:32))-repmat(1:32,[1000000,1]))
