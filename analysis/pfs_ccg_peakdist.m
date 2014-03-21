function [pfsdist,pfsCpow] = pfs_ccg_peakdist(Trial,units,states)

Trial.load('nq');
units = find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19);
nu = numel(units);
states = {'theta','rear&theta','walk&theta'};
numsts = numel(states);

pfs={};
for i = 1:numsts,
    pfs{i}  =  MTAAknnpfs(Trial,units,states{i},0,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end



pfsdist = nan([nu,nu,numsts,numsts]);
pfsCpow = nan([nu,nu,numsts,numsts]);
bns = pfs{1}.adata.binSizes';

for i = 1:numsts
    for j = 1:numsts
        for u = 1:nu
            for o = 1:nu
                pfr = pfs{i}.plot(units(u));
                pfw = pfs{j}.plot(units(o));
                pfr(isnan(pfr)) = 0;
                pfw(isnan(pfw)) = 0;
                C = xcorr2(pfr,pfw);
                cind = [0,0];
                try
                    cind = LocalMinima2(-C,0,bns(1)*2);
                    pfsdist(u,o,i,j) = sqrt(sum((cind-bns).^2));
                    ccind = mat2cell(cind,[1],[1,1]);
                    pfsCpow(u,o,i,j) = C(ccind{:});
                end
            end
        end
    end
end




upfsdist = nan([size(Trial.spk.map,1),size(Trial.spk.map,1),numsts,numsts]);
upfsCpow = nan([size(Trial.spk.map,1),size(Trial.spk.map,1),numsts,numsts]);
upfsdist(units,units,:,:) = pfsdist;
upfsCpow(units,units,:,:) = pfsCpow;

figure,
subplot(121),imagesc(log10(upfsCpow(:,:,3,3)))
subplot(122),imagesc(upfsdist(:,:,3,3))