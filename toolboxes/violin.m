function[h,L,MX,MED,bw]=violin(Y,varargin)
% function[h,L,MX,MED,bw]=violin(Y,varargin)
% 
% violin.m - Simple violin plot using matlab default kernel density estimation
% Last update: 10/2021 
%__________________________________________________________________________
% This function creates violin plots based on kernel density estimation
% using ksdensity with default settings. Please be careful when comparing pdfs
% estimated with different bandwidth!
%
% Differently to other boxplot functions, you may specify the x-position.
% This is usefule when overlaying with other data / plots.
%__________________________________________________________________________
%
% Please cite this function as:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%
%__________________________________________________________________________
%
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
% bw:        Kernel bandwidth. (default []); prescribe if wanted as follows:
%            a) if bw is a single number, bw will be applied to all
%            columns or cells
%            b) if bw is an array of 1xm or mx1, bw(i) will be applied to cell or column (i).
%            c) if bw is empty (default []), the optimal bandwidth for
%            gaussian kernel is used (see Matlab documentation for
%            ksdensity()
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% bw:    bandwidth of kernel
%



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('xlabels'    ,              {{1:numel(Y)}},...
                 'facecolor'  ,              repmat([1, 0.5, 0],[numel(Y),1]),...
                 'linecolor'  ,              'k',...
                 'alpha'      ,              0.5,...
                 'meancolor'  ,              'k',...
                 'mediancolor',              'r',...
                 'bandwidth'  ,              [],...
                 'plotlegend' ,              true,...
                 'plotmean'   ,              true,...
                 'plotmedian' ,              true, ...
                 'xticks'     ,              1:numel(Y) ...
);

[xlabels, facecolor, linecolor, alpha, meancolor, mediancolor, bandwidth, ...
plotlegend, plotmean, plotmedian, xticks] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end

if isempty(find(strcmp(varargin,'bw')))==0
    bandwidth = varargin{find(strcmp(varargin,'bw'))+1}
    if length(bandwidth)==1
        disp(['same bandwidth bw = ',num2str(bandwidth),' used for all cols'])
        bandwidth=repmat(bandwidth,size(Y,2),1);
    elseif length(bandwidth)~=size(Y,2)
        warning('length(bandwidth)~=size(Y,2)')
        error('please provide only one bandwidth or an array of bandwidth with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end


%% Calculate the kernel density
i=1;
for i=1:size(Y,2)
    
    if isempty(bandwidth)==0
        [f, u, bb]=ksdensity(Y{i},'bandwidth',bandwidth(i));
    elseif isempty(bandwidth)
        [f, u, bb]=ksdensity(Y{i});
    end
    
    f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
    
end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------

%% Plot the violins
p = gobjects([1,2]);
i=1;
for i=i:size(Y,2)
    h(i)=fill([F(:,i)+xticks(i);flipud(xticks(i)-F(:,i))],...
              [U(:,i);flipud(U(:,i))],                    ...
              facecolor(i,:),                             ...
              'FaceAlpha',alpha,                          ...
              'EdgeColor',linecolor);


    hold(gca,'on');
    if plotmean == 1
        p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+xticks(i)-i,                 ...
                   interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+xticks(i)-i],...
                  [MX(:,i) MX(:,i)],                                        ...
                  meancolor,                                                ...
                  'LineWidth',1);
    end

    if plotmedian == 1
        p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+xticks(i)-i, ...
                   interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+xticks(i)-i],...
                  [MED(:,i) MED(:,i)],...
                  mediancolor,...
                  'LineWidth',1);
    end
end

%% Add legend if requested
lLabels = {'Mean','Median'};
L=[];
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    L = legend(p(find([plotmean,plotmedian])),lLabels(find([plotmean,plotmedian])))
    set(L,'box','off','FontSize',8)

end

%axis([min(xticks)-0.05*range([min(min(bsxfun(@plus,F,xticks)))), max(xticks)+0.05*range(xticks), min(U(:)) max(U(:))]);
axis('tight');
xlim(xlim+range(xlim)*0.05.*[-1,1]);

%% Set x-labels
set(gca,'XTick',xticks)
box('on');
set(gca,'XTickLabel',xlabels)

%-------------------------------------------------------------------------
