function [uccg,tbins] = unit_ccg(Trial,units,states,spkmode)
%[uccg,tbins] = unit_ccg(Trial,units,states,spkmode)
Trial.load('nq');
[units] = DefaultArgs(varargin,{find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19),'theta'});

nu = numel(units);

Trial.spk.create(Trial,Trial.sampleRate,states,units,'deburst');

units = unique(Trial.spk.clu)';
halfbins = 200;
binsize = 320;
fwin = 31;

[fccg,tbins] = CCG(Trial.spk.res,Trial.spk.clu,binsize,halfbins,Trial.sampleRate,units,'hz');
uccg.ccg = fccg;
tccg = zeros(halfbins*2+1,nu,nu);
cunit = unique(Trial.spk.clu);
cunit = ismember(units,cunit);
tccg(:,cunit,cunit) = fccg;
uccg.fccg = reshape(Filter0(gausswin(fwin)./sum(gausswin(fwin)),tccg),size(tccg));

[uccg.cpow,uccg.tshf] = max(tccg);
uccg.cpow = sq(uccg.cpow);
uccg.tshf = tbins(sq(uccg.tshf));



% $$$ figure
% $$$ subplot(121),imagesc(log10(uccg.cpow)')
% $$$ subplot(122),imagesc(mut')
% $$$ 
% $$$ states = {'theta','rear&theta','walk&theta'};
% $$$ numsts = numel(states);
% $$$ pfs={};
% $$$ for i = 1:numsts,
% $$$     pfs{i}  =  MTAAknnpfs(Trial,[],states{i},0,'numIter',1, ...
% $$$                           'ufrShufBlockSize',0,'binDims',[20, ...
% $$$                         20],'distThreshold',70); ...
% $$$                                                                                                    
% $$$ end




% $$$ unt = [18,91];
% $$$ figure,
% $$$ for s = 1:3,
% $$$ subplot2(2,3,1,s),pfs{s}.plot(pfs{s}.data.clu(unit(1)));
% $$$ subplot2(2,3,2,s),pfs{s}.plot(pfs{s}.data.clu(unit(2)));
% $$$ end
% $$$ uic = mat2cell(unit,[1],[1,1]);
% $$$ muc(uic{:})
% $$$ mut(uic{:})
% $$$ 
% $$$ xyc = round(get(gca,'CurrentPoint'));
% $$$ unit = [18,91];
% $$$ figure,
% $$$ for s = 1:3,
% $$$ subplot2(2,3,1,s),pfs{s}.plot(88);
% $$$ subplot2(2,3,2,s),pfs{s}.plot(35);
% $$$ end
% $$$ uic = mat2cell(unit,[1],[1,1]);
% $$$ muc(uic{:})
% $$$ mut(uic{:})
% $$$ 
% $$$ 
% $$$ figure
% $$$ nu = 10;
% $$$ unit = 22:31;
% $$$ for u = 1:nu,
% $$$ for o = 1:nu,
% $$$ subplot2(nu,nu,u,o);
% $$$ bar(t,uccg(:,unit(u),unit(o)));
% $$$ end
% $$$ end