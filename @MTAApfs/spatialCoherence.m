function rho = spatialCoherence(Pfs,units,maskFlag)
% rho = spatialCoherence(Pfs,units);
% The correlation between a list of firing rates in each
% pixel and a corresponding list of firing rates averaged
% over the 8 nearestneighbors of each pixel. 
% muller_kubie_1989
if isempty(units),
    units = Pfs.data.clu;
end
rho = nan(numel(units),1);
for u = units(:)',
    rm = Pfs.plot(u,'mazeMaskFlag',maskFlag);
    cmat = [1,1,1;1,0,1;1,1,1];
    c = conv2(rm,cmat,'valid')./8;
    rmss = rm(2:end-1,2:end-1);
    ind = ~isnan(c)&~isnan(rmss);
    rho(u==units) = corr(c(ind),rmss(ind));
end

