function unit_dist(Trial,varargin)
Trial.load('nq');
[units,states] = DefaultArgs(varargin,{find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19),{'theta','rear&theta','walk&theta'}});
numsts = numel(states);
nu = numel(units);

pfs={};
for i = 1:numsts,
    pfs{i}  =  MTAAknnpfs(Trial,units,states{i},0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
end

[uccg,tbins] = unit_ccg(Trial,units,states);
[pccg] = pfs_ccg_peakdist(Trial,units,states);


unt = 38;
figure(1293)
clf,
udtr = sqrt(pccg.dist(:,unt,1,2).^2+(1./uccg.cpow(:,unt)).^2);
udtw = sqrt(pccg.dist(:,unt,1,3).^2+(1./uccg.cpow(:,unt)).^2);
udt = sqrt(pccg.dist(:,unt,1,1).^2+(1./uccg.cpow(:,unt)).^2);
udr = sqrt(pccg.dist(:,unt,2,2).^2+(1./uccg.cpow(:,unt)).^2);
udw = sqrt(pccg.dist(:,unt,3,3).^2+(1./uccg.cpow(:,unt)).^2);
subplot(331),pfs{1}.plot(units(unt));
subplot(332),pfs{2}.plot(units(unt));
subplot(333),pfs{3}.plot(units(unt));
subplot(334),plot(pccg.dist(:,unt,1,1),1./uccg.cpow(:,unt),'.'),
subplot(335),plot(pccg.dist(:,unt,2,2),1./uccg.cpow(:,unt),'.'),
subplot(336),plot(pccg.dist(:,unt,3,3),1./uccg.cpow(:,unt),'.'),
subplot(337),hist(udt,30)
subplot(338),hist(udr,30)
subplot(339),hist(udw,30)

figure(3884),clf
plot(udtw,udtr,'.')
hold on,line([0,45],[0,45])
%subplot(324),plot(pccg.dist(:,unt,2,3),1./uccg.cpow(:,unt),'.'),


figure(324244),clf,plot(pccg.dist(:,unt,3,3),pccg.dist(:,unt,2,2),'.')
hold on,line([0,45],[0,45])

figure(324244),clf,unt = 55;
subplot(121),plot(pccg.dist(:,unt,2,2),log10(pccg.cpow(:,unt,2,2)),'.')
subplot(122),plot(pccg.dist(:,unt,3,3),log10(pccg.cpow(:,unt,3,3)),'.')


