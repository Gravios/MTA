function cpnts = ClusterH2(Data,Xbins,Ybins)
%function cpnts = ClusterH2(Data,Xbins,Ybins)
%
% See hist2
%
dind = nniz(Data);
hfig = figure;
[~,xb,yb,pos] = hist2(Data(dind,:),Xbins,Ybins);
nind = nniz(pos);
plot(xb(pos(nind,1)),yb(pos(nind,2)),'.');
hold on,
hist2(Data,Xbins,Ybins);
nind(nind==1) = ClusterPP(hfig);
cpnts = zeros([size(Data,1),1]);
cpnts(dind) = nind;