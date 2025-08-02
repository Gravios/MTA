


Trial = MTATrial.validate('jg05-20120312.cof.all');
[units,nq] = select_units(Trial,'int');


sampleRate = 250;
xyz = preproc_xyz(Trial,'trb',sampleRate);


hvang = filter(copy(xyz),'ButFilter',3,1.5,'low');
xycoor = cat(2,...
             hvang(:,'spine_upper',[1,2])-hvang(:,'bcom',[1,2]),...
             hvang(:,'nose',[1,2])-hvang(:,'hcom',[1,2]));
hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% Positive: CCW (Left)     Negative: CW (Right)
hvang.data = circ_dist(circshift(circ_dist(hvang.data(:,2),hvang.data(:,1)),-10),...
                          circshift(circ_dist(hvang.data(:,2),hvang.data(:,1)),+10));

hbang = filter(copy(xyz),'ButFilter',3,30,'low');    
xycoor = cat(2,...
             hbang(:,'spine_upper',[1,2])-hbang(:,'bcom',[1,2]),...
             hbang(:,'nose',[1,2])-hbang(:,'hcom',[1,2]));
hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
hbang.data = circ_dist(hbang.data(:,2),hbang.data(:,1));




phz = phase( load( copy( Trial.lfp ), Trial, 69), [6, 12]);

phz.data = unwrap(phz.data);
resample(phz,xyz);
phz.data = mod(phz.data+pi,2*pi)-pi;

hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(hbang.data-0.2),hbangBinEdges);

phzBinEdges = linspace(-pi,pi,12);
phzBinCenters = mean([phzBinEdges(2:end); phzBinEdges(1:end-1)]);
edy = numel(phzBinCenters);
phzBinInd = discretize(phz.data,phzBinEdges);

posxBinEdges = linspace(-400,400,10);
posxBinCenters = mean([posxBinEdges(2:end); posxBinEdges(1:end-1)]);
x = numel(posxBinCenters);
posxBinInd = discretize(xyz(:,'hcom',1),posxBinEdges);

posyBinEdges = linspace(-400,400,10);
posyBinCenters = mean([posyBinEdges(2:end); posyBinEdges(1:end-1)]);
y = numel(posyBinCenters);
posyBinInd = discretize(xyz(:,'hcom',2),posyBinEdges);




states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
stcm = stc2mat(Trial.stc,xyz,states);

spk = Trial.spk.copy;
spk.create(Trial,xyz.sampleRate,'theta',units,'none');

ind = stcm(:,1) == 1 & any(stcm(:,[3,5,6]),2);
ind = ind & nniz([hbangBinInd,phzBinInd,posxBinInd,posyBinInd]);
occ = accumarray([hbangBinInd(ind),phzBinInd(ind),posxBinInd(ind),posyBinInd(ind)],1,[edx,edy,x,y],@sum)./sampleRate;

figure
sax = reshape(tight_subplot(y,x,0.001,0,0),[x,y])';
for u = 1:numel(units),
res = spk(units(u));
res = res(nniz([phzBinInd(res),hbangBinInd(res),ind(res),posxBinInd(res),posyBinInd(res)]));
out = accumarray([hbangBinInd(res),phzBinInd(res),posxBinInd(res),posyBinInd(res)],1,[edx,edy,x,y],@sum);
for xi = 1:x,
    for yi = 1:y,
        axes(sax(yi,xi));
        imagesc(out(:,:,xi,yi)'./occ(:,:,xi,yi)');
        
    end
end
af(@(s) caxis(s,[0,prctile(out(:),[99.99])]), sax)
axes(sax(1));
title(num2str(units(u)));
colorbar();
waitforbuttonpress();
end