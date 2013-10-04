function batchFigFontSizeChangingFunction(varargin)
[fig_path,fontSize] = DefaultArgs(varargin,{'/data/homes/gravio/figures/asFIG-20130512/',20});

files = dir(fig_path);
re = ['\.fig'];
figFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
for fig = figFileList,
    hfig = open([fig_path fig{1}]);
    children = get(hfig,'Children')';
    for child = children,
        set(child,'FontSize',fontSize);
    end
    print(gcf,'-dpsc2',[fig_path fig{1}(1:end-3) 'eps']);
end
