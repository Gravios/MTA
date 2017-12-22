function plot_features_with_stc(features,stc,varargin)
defargs = struct('sampleRate',   1,       ...
                 'plotFunction', 'plot',  ...
                 'lims',         [],      ...
                 'fetLabels',    {{}}     ...
);
[sampleRate,plotFunction,lims,fetLabels] = DefaultArgs(varargin,defargs,'--struct');


figure();

% DISPLAY features
subplot2(3,1,[1,2],1);
ts = [1:size(features,1)]./features.sampleRate;
switch plotFunction
  case 'plot'
    plot(ts,features(:,:));
    if ~isempty(lims),  ylim(lims);  end
    if ~isempty(fetLabels),legend(fetLabels{:});end
    
  case 'imagesc'
    imagesc(ts,1:size(features,2),features(:,:)');axis xy
    if ~isempty(lims),  caxis(lims);  end
  
  case 'imagescnan'
    
    
end    

% DISPLAY behavioral states
subplot2(3,1,3,1);
plotSTC(stc,sampleRate);

% LINK the x axes of the two supblots 
linkaxes(findobj(gcf,'Type','axes'),'x')


