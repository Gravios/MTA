function ppf(bin1,bin2,rateMap)
imagesc(bin1,bin2,rateMap)
colorbar
colormap('default')
cc = colormap;
cc(1,:) = [0 0 0];
colormap(cc)
imagesc(bin1,bin2,rateMap')
if ~isempty(rateMap)&~isempty(bin1)&~isempty(bin2),
text(bin1(1)+10,bin2(end)-10,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
end
axis xy
