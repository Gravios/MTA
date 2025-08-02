function plot_stcs(varargin)

nsts = numel(varargin);
figure,
for i = 1:nsts,
    subplot(nsts,1,i);
    plotSTC(varargin{i});
end
linkaxes(findobj(gcf,'Type','axes'),'x');