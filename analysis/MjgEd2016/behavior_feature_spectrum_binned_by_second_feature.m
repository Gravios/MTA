function behavior_feature_spectrum_binned_by_second_feature(Trial,varargin)
%behavior_feature_spectrum_binned_by_second_feature

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(...
    'Feature',        [],...
    'frequencyBins',    [],...
    'mode',           {{'bspeed'}},...
    'state',          'walk',...
    'stcMode',        'NN0317',...
    'marker',         'spine_lower',...
    'figureHandleNum',gen_figure_id);
%--------------------------------------------------------------------------------------------------

[Feature,frequencyBins,mode,state,stcMode,marker,figureHandleNum] = DefaultArgs(varargin,defargs,'--struct');

if isempty(Feature),
    [Feature,frequencyBins] = fet_rhm(Trial,[],'mtchglong');
    Feature.data  = log10(Feature.data);
    Feature.data(Feature<-9) = nan;
    Feature.data(nniz(Feature.data))=nan;
elseif ischar(Feature),
    
end


% MAIN ---------------------------------------------------------------------------------------------
if ischar(mode), mode = {mode}; end

Trial = MTATrial.validate(Trial);

Stc = Trial.load('stc',stcMode);
state = Stc{state};

% Compute the duration of each walk period
wdur = diff(state.data,1,2);
% remove periods which are shorter than 1 second
state.data(wdur<90,:)=[];


xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

% Figure settings
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8);
hfig = figure(figureHandleNum); clf;
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,(4+2)*numel(mode)+2,2*(3+2)+2]);
set(hfig,'PaperPositionMode','auto');


vel = xyz.vel(marker,[1,2]);
vel.resample(Feature);
vnn = nniz(vel);

% $$$ Feature.data = (Feature.data-repmat(nanmean(Feature(vnn,:)),[Feature.size(1),1]))...
% $$$            ./repmat(nanstd(Feature(vnn,:)),[Feature.size(1),1]);


for m = 1:numel(mode),

    switch mode{m}
      case 'height'
        vel = MTADxyz('data',xyz(:,'head_front',3),'sampleRate',xyz.sampleRate);
        vel.resample(Feature);
        vel.data = log10(vel.data);
        vnn = nniz(vel);
        vhlim = [0, 2.5];
        vh_label = 'Height(H)';
        vh_units = 'log10(mm)';
      
      case 'hangle'
        ang = create(MTADang,Trial,xyz);
        vel = MTADxyz('data',ang(:,'head_back','head_front',2),'sampleRate',xyz.sampleRate);
        vel.resample(Feature);
        vnn = nniz(vel);
        vh_label = 'Pitch(H)';
        vh_units = 'rad';
        vhlim = [-pi/2, pi/2];        
      
      case 'bspeed'
        vel = xyz.vel(marker,[1,2]);
        vel.resample(Feature);
        %vel.data = log10(vel.data);
        vnn = nniz(vel);
        vhlim = prctile(vel(state),0:5:100);
        %vhlim = [0,2.2];
        vh_label = 'Speed(B)';
        vh_units = 'cm/s';
    
    end


    sFeature = Feature(state,:);
    svel = vel(state,:);

    vbins = 10;
    %vedgs = linspace(vhlim(1),vhlim(2),vbins);
    vbins = numel(vhlim);
    vedgs = vhlim;
    [N,vbs] = histc(svel,vedgs);


    mrv = nan(numel(N),Feature.size(2));
    for f =1:Feature.size(2),
        mrv(:,f) = accumarray(vbs(nniz(vbs)),sFeature(nniz(vbs),f),[vbins,1],@nanmean);
    end
    mrv(N<10,:) = nan;


    axes('Units','centimeters',...
         'Position',[(4+2)*(m-1)+2,(3+2)*2-(4)+1,4,3])

    
    imhand = imagesc(vedgs,frequencyBins,bsxfun(@plus,mrv,abs(min(mrv(:))))');
    axis('xy');
    set(imhand,'tag', [mfilename,'-',mode{m}])
    title([Feature.label ' mean PSD given ' vh_label])
    xlabel([ vh_label ' (' vh_units ')']);
    ylabel('Frequency Hz')
    colormap jet
    caxis(prctile(reshape(get(get(gca,'Children'),'CData'),[],1),[2,98]))
    h=colorbar('EastOutside');
    apos = get(gca,'Position');
    set(h,'Units','centimeters');
    set(h,'Position',[apos(1)+4.1,apos(2),0.3,3])


    axes('Units','centimeters',...
                 'Position',[(4+2)*(m-1)+2,(3+2)*1-(4)+1,4,3])
            %subplot2(numel(mode),2,m,2);
    bar(vedgs,N,'histc')
    title(['Marginal Distrb of ' mode{m}]);
    xlabel([ vh_label ' (' vh_units ')']);
    ylabel('count');
    xlim(vedgs([1,end]));

end


tbox = uicontrol('Style','text',...
          'Parent',hfig,...
          'String',[Trial.filebase ' : ' state.label],...
          'Units','centimeters',...
          'Position',[2,(3+2)*2+1,4*numel(mode),.5])

tbox.Position([3,4]) = tbox.Extent([3,4]);

figPath = fullfile(Trial.spath,'figures');
if ~exist(figPath), mkdir(figPath),end
savefig(hfig,fullfile(figPath,[Feature.label,'_psd_distrib_',strjoin(mode,'_'),'.fig']));
print(gcf,'-depsc2',fullfile(figPath,[Feature.label,'_psd_distrib_',strjoin(mode,'_'),'.eps']));


%---------------------------------------------------------------------------------------------------
