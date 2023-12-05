function render_subplot_by_name(hfig,document,subplotPath,subplotName)

% partsPath
% document
% hfig
% subplotName

delete(findobj(hfig,'Tag',subplotName));

widthScale  = document.body.width;
heightScale = document.body.height;

tempfig = openfig(fullfile(subplotPath,[subplotName,'.fig']));
tempfig.Units='normalized';
ax = copyobj(tempfig.Children,hfig);

% GET subplot handle
sax = findobj(ax,'Tag',subplotName)
% REPOSITION subplot
for a = 1:numel(sax)
    sax(a).Position = [document.subplots.(subplotName).x/widthScale, ...
                    1 - document.subplots.(subplotName).y/heightScale - document.subplots.(subplotName).height/heightScale,...
                    document.subplots.(subplotName).width/widthScale,...
                    document.subplots.(subplotName).height/heightScale];
end
close(tempfig);

disp(subplotName);
