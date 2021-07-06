function report_unit_summary(Trial,varargin)
% function report_placefield_summary
% 
% description: plot a variety of measures for units with place fields 
%
% Mod:20171211: reorganize the signature of multiple functions to include units
%               as the second input option.
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('figDir',                '/storage/gravio/figures/analysis/units',        ...
                 'pitchReferenceTrial',   'Ed05-20140529.ont.all',                               ...
                 'marker',                'hcom',                                                ...
                 'thetaRef',              1,                                                     ...
                 'overwrite',             false                                                  ...
);
[figDir,pitchReferenceTrial,marker,thetaRef,overwrite] = DefaultArgs(varargin,defargs,'--struct');

configure_default_args();


% $$$ Trial = MTATrial.validate('FS03-20201222.cof.all');
% $$$ figDir = '/storage/gravio/figures/analysis/units';
create_directory(figDir);
create_directory(fullfile(figDir,Trial.filebase));
units = Trial.spk.map(:,1)';
stc = Trial.load('stc','msnn_ppsvd_raux');



states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
statesLabels = {'theta','rear','hloc','hpause','lloc','lpause'};
numStates = numel(states);

Trial.load('nq');

pft = pfs_2d_theta(Trial,units,'overwrite',true);
pfs = pfs_2d_states(Trial,units,'msnn_ppsvd_raux',states,'overwrite',true);
bfs = compute_bhv_ratemaps(Trial,units,'overwrite',true);

[accg,tbins] = autoccg(Trial);


tppSampleRate = 250;
spkw = Trial.spk.copy();
spkw.load_spk(Trial,tppSampleRate);

xyz = preproc_xyz(Trial,'trb');
xyz.resample(tppSampleRate);

% CREATE lowpass filtered xyz object
% COMPUTE basis vector aligned to the head
fxyz = filter(copy(xyz),'ButFilter',3,20,'low');    
hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);


%drz = compute_drz(Trial,units,pfs{1},[],[],'hcom',[],xyz);
%hrz = compute_hrz(Trial,units,pfs{1},[],[],'hcom',[],xyz);
ghz = compute_ghz(Trial,units,pfs{1},[],[],'hcom',[],xyz);    
gdz = compute_gdz(Trial,units,pfs{1},[],[],'hcom',[],xyz);    
ddz = compute_ddz(Trial,units,pfs{1},[],[],'hcom',[],xyz);    

phz = load_theta_phase(Trial,xyz,2,0);

global MTA_PROJECT_PATH
% MAKE xyhb ratemap mask
% $$$ mask = repmat(maskXY,[1,1,size(maskbhv)]).*repmat(permute(maskbhv,[3,4,1,2]),[size(maskXY),1,1]);
%save(fullfile(MTA_PROJECT_PATH,'analysis','pfsHB_mask.mat'),'mask','-v7.3');
% LOAD xyhb ratemap mask
ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsHB_mask.mat'));
maskBhv = ds.mask;
maskXY = create_tensor_mask(pfs{1}.adata.bins);


% TODO generate saved mask 

axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','auto'};
nanColor = [0.3,0.3,0.3];

[hfig,fig,fax,sax] = set_figure_layout(figure(666000),'A4','landscape',[],2.5,2.5,0.2,0.2);
%for unit = units

unit = 1;
while unit~=-1    
    clf();
disp(unit);    
sax = gobjects([0,1]);


[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, -1);
FigInfo = uicontrol('Parent',hfig,                                                        ...
                    'Style','text',                                                       ...
                    'String',{['Unit: ',num2str(unit)],                                   ...
                              Trial.filebase,                                             ...
                              ['stcMode: ',   Trial.stc.mode],                            ...
                              ['eDist:   ',   num2str(Trial.nq.eDist(unit))],             ...
                              ['Refrac:  ',   num2str(log10(Trial.nq.Refrac(unit)))],     ...
                              ['SNR:     ',   num2str(Trial.nq.SNR(unit))],               ...
                              ['AmpSym:  ',   num2str(Trial.nq.AmpSym(unit))],            ...
                              ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]          ...
                             },                                                           ...
                    'Units','centimeters',                                                ... units
                    'Position',[fig.page.xpos(xind)+xOffSet,                              ... x
                                fig.page.ypos(yind)+yOffSet,                              ... y
                                6,                                                        ... w
                                3.5]                                                      ... h
);

% SECTION place field and bhvfield


mrate = max(cell2mat(cf(@(p) p.maxRate(unit),pfs)));

% PLOT bhvfield
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 3, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height]);
hold(sax(end),'on');
plot(bfs,unit,'mean','text',[0,mrate],false,[],false,[],@jet,maskBhv,nanColor);
tax = findobj(gca(),'Type','Text');
tax.Position = [-0.15,1.5,0];
sax(end).YTickLabel = {};
sax(end).XTickLabel = {};
xlim([-1.65,0.5]);
ylim([-0.5,1.65]);


% PLOT placefields by state
for sts = 1:numStates
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, sts+3, 0);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    hold(sax(end),'on');
    plot(pfs{sts},unit,'mean','text',[0,mrate],true,[],false,[],@jet,maskXY,nanColor);
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    title(statesLabels{sts});
end



% PLOT accg    
[yind, yOffSet, xind, xOffSet] = deal(2, -0.5, 1, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*1.5,                        ...
                              fig.subplot.height]);
bar(tbins,accg(:,unit));axis tight;
sax(end).YTickLabel = {};
sax(end).XTickLabel = {};


% PLOT theta phase distribution    
[yind, yOffSet, xind, xOffSet] = deal(2, -0.5, 2, fig.subplot.width-0.5);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height]);
res = spkw.res(spkw.clu==unit);
tper = [stc{'t',tppSampleRate}];
res = res(WithinRanges(res,tper.data));
rose(phz(res,1));axis tight;
sax(end).YTickLabel = {};
sax(end).XTickLabel = {};


% PLOT spk waveform
[yind, yOffSet, xind, xOffSet] = deal(5, 0, 1, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*1.5,                    ...
                              fig.subplot.height*3]);
uResw = spkw.res(spkw.clu==unit);
uSpkw = spkw.spk(spkw.clu==unit,:,:);
[~,sInd] = SelectPeriods(uResw,[stc{states{1},tppSampleRate}],'d',1,0);
mspk = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
sspk = sq(std(uSpkw(sInd,:,:)))';    

hold('on');
plot(mspk,'b');
plot(mspk+sspk,'r');    
plot(mspk-sspk,'r');        
xlim([0,size(spkw.spk,3)]);
sax(end).YTickLabel = {};
sax(end).XTickLabel = {};
box(sax(end),'on');





% PLOT phase precession
for sts = 1:numStates,
    % GET current state
    % SELECT spikes within current state 
    % SKIP if spike count is less than 50
    state = [stc{states{sts},spkw.sampleRate}];
    res = spkw(unit);
    res = res(WithinRanges(res,state.data));
    res(abs(ddz(res,unit))>350) = [];
    if numel(res) >10,
        res(res>xyz.size(1))=[];
        %drzspk = drz(res,u);
        %hrzspk = hrz(res,unit);                
        ghzspk = ghz(res,unit);
        gdzspk = gdz(res,unit);
        %ddzspk = ddz(res,u);
        phzspk = phz(res,1);
        %gind = ~isnan(drzspk)&~isnan(phzspk)&~isnan(ghzspk)&~isnan(gdzspk);
        %gind = ~isnan(hrzspk)&~isnan(phzspk)&~isnan(ghzspk);
        gind = ~isnan(gdzspk)&~isnan(phzspk)&~isnan(ghzspk);
    else
        res = [];
        drzspk=[];
        ddzspk=[];
        phzspk=[];
        ghzspk=[];
        gdzspk=[];
        gind=[];
    end
    [yind, yOffSet, xind, xOffSet] = deal(2, 0, sts+3, 0);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                        fig.page.ypos(yind)+yOffSet,              ...
                        fig.subplot.width,                        ...
                        fig.subplot.height]);
    hold(sax(end),'on');
    if sum(gind)>10,
        plot([ghzspk(gind);ghzspk(gind)],...
             [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
             '.','MarkerSize',4);
        xlim([-1,1]);
        ylim([0,720]);        
    end        
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};

    
    [yind, yOffSet, xind, xOffSet] = deal(3, 0, sts+3, 0);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                        fig.page.ypos(yind)+yOffSet,              ...
                        fig.subplot.width,                        ...
                        fig.subplot.height]);
    hold(sax(end),'on');
    if sum(gind)>10,
        plot([gdzspk(gind);gdzspk(gind)],...
             [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
             '.','MarkerSize',4);
        xlim([-1,1]);
        ylim([0,720]);
    end
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    

        % GET placefield center
        % COMPUTE position of the placefield center relative to head basis
        [mxr,mxp] = pfs{1}.maxRate(unit);
        pfsCenterHR = MTADfet.encapsulate(Trial,                                               ...
                                          multiprod(bsxfun(@minus,mxp,sq(xyz(:,'nose',[1,2]))),...
                                                    hvec,2,[2,3]),                             ...
                                          tppSampleRate,                                       ...
                                          'placefield_center_referenced_to_head',              ...
                                          'pfsCenterHR',                                       ...
                                          'p'                                                  ...
                                          );
        
    
        %res = res(WithinRanges(res,get([stc{'x+p&t'}],'data'))&cind(res));
        
        % MEAN spike theta phase at location of place field center relative to head (Wait ... what?)
    [yind, yOffSet, xind, xOffSet] = deal(4, 0, sts+3, 0);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                        fig.page.ypos(yind)+yOffSet,              ...
                        fig.subplot.width,                        ...
                        fig.subplot.height]);
    hold(sax(end),'on');
    if sum(gind)>10&sum(sqrt(sum(pfsCenterHR(res,:).^2,2))<350)>10,
        myPhz = phz(res,1);
        myPhz(myPhz>pi) = myPhz(myPhz>pi)-2*pi;
        scatter(pfsCenterHR(res,1),pfsCenterHR(res,2),5,...
                myPhz,'filled');
        colormap(sax(end),'hsv');
        caxis([-pi,pi]);
        set(gca,'Color',[0.75,0.75,0.75]);
        ylim([-300,300]);
        xlim([-300,300]);
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};
        
    end
    if sts==1,
        ylabel('egocentric');
    end
    

        % MEAN spike theta phase at location of place field center relative to head (Wait ... what?)
    [yind, yOffSet, xind, xOffSet] = deal(5, 0, sts+3, 0);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                        fig.page.ypos(yind)+yOffSet,              ...
                        fig.subplot.width,                        ...
                        fig.subplot.height]);
    hold(sax(end),'on');
        

        if sum(gind)>10&sum(sqrt(sum(pfsCenterHR(res,:).^2,2))<350)>10,
            bins = linspace(-300,300,15);
            hind = discretize(pfsCenterHR(res,:),bins);
            mtph = accumarray(hind(nniz(hind),:),                                ... subs
            phz(res(nniz(hind)),1),... vals
            repmat(numel(bins),[1,2]),                         ... size
            @circ_mean);                                         % func
            mtph(mtph==0) = nan;
            pax = pcolor(bins,bins,mtph');
            pax.EdgeColor = 'none';
            caxis(sax(end),[-pi,pi]);
            colormap(sax(end),'hsv');
            if sts == numStates
                cax = colorbar(sax(end));
                cax.Units = 'centimeters';
                cax.Position(1) = sum(sax(end).Position([1,3])) + 0.1;
            end
            %imagescnan({bins,bins,mtph'},[0,2*pi],'circular',true,'colorMap',@hsv);
            xlim([-300,300])
            ylim([-300,300])            
            axis('xy');
            sax(end).YTickLabel = {};
            sax(end).XTickLabel = {};

        end

    if sts==1,
        ylabel('egocentric');
    end
    

end
    pause(0.01);
    
    FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit,'%04.f')];
    %print(hfig,'-depsc2',fullfile(figDir,Trial.filebase,[FigName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,Trial.filebase,[FigName,'.png']));

unit = figure_controls(hfig,unit,units(:)',true);

end
