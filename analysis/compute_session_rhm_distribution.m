function compute_session_rhm_distribution(Trial,varargin)

% DEFARGS ----------------------------------------------------------------------
defargs = struct('mode',    {{'height','hangle'}},                           ...
                 'state',   'loc',                                           ...
                 'stcMode', 'NN0317',                                        ...
                 'marker',  'spine_lower',                                   ...
                 'figHnum', 20160926                                         ...
);

[mode,state,stcMode,marker,figHnum] = DefaultArgs(varargin,defargs,'--struct');
if ischar(mode), mode = {mode}; end
% END DEFARGS ------------------------------------------------------------------



% DEFVARS ----------------------------------------------------------------------
Trial = MTATrial.validate(Trial);
Stc = Trial.load('stc',stcMode);
state = Stc{state};
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
% END DEFVARS ----------------------------------------------------------------------


    
% MAIN -------------------------------------------------------------------------    

% figure options
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8);
hfig = figure(figHnum);
clf;
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,(4+2)*numel(mode)+2,2*(3+2)+2]);
set(hfig,'PaperPositionMode','auto');

FigPath = fullfile(Trial.spath,'figures');
if ~exist(FigPath), mkdir(FigPath),end


% Get the rhm feature
[rhm,fs] = fet_rhm(Trial,[],'mtchglong');


% log10 and unity
rhm.data  = log10(rhm.data);
rhm.data(rhm<-9) = nan;
rhm.data(nniz(rhm.data))=nan;

vel = xyz.vel(marker,[1,2]);
vel.resample(rhm);
vnn = nniz(vel);

rhm.data = (rhm.data-repmat(nanmean(rhm(vnn,:)),[rhm.size(1),1]))...
           ./repmat(nanstd(rhm(vnn,:)),[rhm.size(1),1]);


% loop through modes 
for m = 1:numel(mode),

    switch mode{m}
      case 'height'
        vel = MTADxyz('data',xyz(:,'head_front',3),'sampleRate',xyz.sampleRate);
        vel.resample(rhm);
        vel.data = log10(vel.data);
        vnn = nniz(vel);
        vhlim = [0, 2.5];
        vh_label = 'Height(H)';
        vh_units = 'log10(mm)';
      
      case 'hangle'
        ang = create(MTADang,Trial,xyz);
        vel = MTADxyz('data',ang(:,'head_back','head_front',2),'sampleRate',xyz.sampleRate);
        vel.resample(rhm);
        vnn = nniz(vel);
        vh_label = 'Pitch(H)';
        vh_units = 'rad';
        vhlim = [-pi/2, pi/2];        
      
      case 'bspeed'
        vel = xyz.vel(marker,[1,2]);
        vel.resample(rhm);
        vel.data = log10(vel.data);
        vnn = nniz(vel);
        vhlim = [-3,2.2];
        vh_label = 'Speed(B)';
        vh_units = 'log10(cm/s)';
    
    end


    srhm = rhm(state,:);
    svel = vel(state,:);

    vbins = 100;
    vedgs = linspace(vhlim(1),vhlim(2),vbins);
    [N,vbs] = histc(svel,vedgs);


    mrv = nan(numel(N),rhm.size(2));
    for f =1:rhm.size(2),
        mrv(:,f) = accumarray(vbs(nniz(vbs)),srhm(nniz(vbs),f),[vbins,1],@nanmean);
    end
    mrv(N<20,:) = nan;


    axes('Units','centimeters',...
         'Position',[(4+2)*(m-1)+2,(3+2)*2-(4)+1,4,3])

    
    imhand = imagesc(vedgs,fs,bsxfun(@plus,mrv,abs(min(mrv(:))))');
    axis('xy');
    set(imhand,'tag', [mfilename,'-',mode{m}])
    title(['RHM mean PSD given ' vh_label])
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

end


tbox = uicontrol('Style','text',...
          'Parent',hfig,...
          'String',[Trial.filebase ' : ' state.label],...
          'Units','centimeters',...
          'Position',[2,(3+2)*2+1,4*numel(mode),.5])
tbox.Position([3,4]) = tbox.Extent([3,4]);


savefig(hfig,fullfile(FigPath,['RHM_psd_distrib_',strjoin(mode,'_'),'.fig']));
print(gcf,'-depsc2',fullfile(FigPath,['RHM_psd_distrib_',strjoin(mode,'_'),'.eps']));

% END MAIN -------------------------------------------------------------------------    
