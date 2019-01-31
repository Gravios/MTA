% req20190122
%    Tags: behavioral information time shifted
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: theta modulation of bhv rate maps
%    Description: Compilation of placefield information and bhv theta decomposition stuff
%    Protocol:
%    Figures: /storage/gravio/figures/analysis/placefields_bhv

global ANAL_PARM

% FIGURE Opts ----------------------------------------------------------------
FigDir = create_directory('/storage/gravio/figures/analysis/placefields_bhv');

pageHeight = 21.0;
pageWidth  = 29.7;

pwidth = 1.5;
pheight = 1.5;

xpad = 0;
ypad = 0;

xpos = 3.5:(pwidth+xpad):pageWidth;
ypos = fliplr(0.5:(pheight+xpad):pageHeight-3.7);

hfig = figure(666002);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';
hfig.PaperType = 'A4';
hfig.PaperUnits = 'centimeters';

% DATA Load ----------------------------------------------------------------
MjgER2016_load_data();
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%

% ANALYSIS Opts
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

pft = cf(@(Trial)  pfs_2d_theta(Trial),  Trials);
pfs = cf(@(Trial) MTAApfs(Trial,'tag','hbpptbpFS1v3'), Trials);
pfb = cf(@(Trial,unitSubset)  pfs_2d_states(Trial),  Trials,units);
pftp = cf(@(Trial) MTAApfs(Trial,'tag','tp'), Trials);

if ~exist('pfd','var'), ...
        [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
        req20180123_ver5(Trials,[],'13');  
end

req20181119();% -> eigVals, eigVec, eigScrs, eigVars, eSpi, FSrBhvThp
req20190122();% -> bhvInfoTS, tshifts

nanColor = [0.1,0.1,0.1];


% FIGURE : Behavioral theta phase decomposition Examples -------------------------------------------%

for t = 20:23;
    Trial = Trials{t}; 
    unitSubset = units{t};
    Trial.load('nq');
    
    % TEMPORARY
    pftp = cell([1,23]);
    pftp{t} = MTAApfs(Trials{t},'tag','tp');

    create_directory(fullfile(FigDir,Trial.filebase));
    disp(['Processing Trial: ' Trial.filebase]);







    % AXES SAMPLE
% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
% $$$                  'FontSize', 8);

    phzOrder = [10,12,14,16,2,4,6,8];
    phzLables = {'','90','','180','','270','',''};



    for unit = unitSubset,
        clf();        

        usi = ismember(cluSessionMap,[t,unit],'rows');
        sp = gobjects([1,0]);        
        fax = axes('Position',  [0,0,1,1],      ...
                   'Visible',   'off',          ...
                   'Units',     'centimeters',  ...
                   'FontSize',  8,              ...
                   'LineWidth', 1);
        xlim(fax,[0,hfig.Position(3)]);
        ylim(fax,[0,hfig.Position(4)]);
        
        % GET space and theta phase jointly restricted bhv rate maps
        rmap = plot(pfs{t},unit,1,'colorbar',[],false);
        rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs{t}.adata.binSizes(2);
        zmap = reshape(permute(rmap,[1,3,2]),[],pfs{t}.adata.binSizes(2));
        zmap(repmat(~rmapNind(:),[1,pfs{t}.adata.binSizes(2)])) = 0;            
        zmap(zmap(:)<0) = 0;            
        rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs{t}.adata.binSizes(2))','omitnan')~=0,...
                                    size(rmap,1),size(rmap,3));      
        mrate= prctile(rmap(:),99.9);    

        if isnan(mrate),continue,end



        % PLOT theta restricted spatial rate map
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(1),ypos(1),pwidth,pheight],...
                         'FontSize', 8);
        plot(pft{t},unit,'mean','text',mrate,'colorMap',@jet,'nanColor',nanColor);
        ylabel(['Unit: ',num2str(unit)]);
        set(sp(end),'YTickLabels',{});
        set(sp(end),'XTickLabels',{});    
        title({'Spatial','Rate Map'});
        % REFERENCE bar physical space
        axes(fax);
        line([sp(end).Position(1),sum(sp(end).Position([1,3]).*[1,0.5])],...
             sp(end).Position(2).*[1,1]-0.1,'LineWidth',1);
        text(sp(end).Position(1),sp(end).Position(2)-0.3,'50 cm','FontSize',8);

        drawnow();

        % PLOT theta and space restricted behavior rate map
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(1),ypos(3),pwidth,pheight],...
                         'FontSize', 8);
        plot(pfd{t},unit,'mean','text',mrate,false,'colorMap',@jet,'nanColor',nanColor);
        set(sp(end),'YTickLabels',{});
        set(sp(end),'XTickLabels',{});    
        ylabel('Body Pitch');
        xlabel('Head Pitch');
        title({'Behavioral','Rate Map'});
        % REFERENCE bar behavioral space
        axes(fax);
        line(sum(sp(end).Position([1,3])).*[1,1]+0.1,...
             [sp(end).Position(2),sum(sp(end).Position([2,4]).*[1,11/28])],'LineWidth',1);
        text(sum(sp(end).Position([1,3]))+0.3,sp(end).Position(2),'1 rad','FontSize',8,'Rotation',90);





        % PLOT behaviorally restricted spatial rate maps
        stateOrder = [4,3,7,2,6];
        for s = 1:5
            sp(end+1) = axes('Units','centimeters',...
                             'Position',[xpos(2+s)+0.2*(s-1),ypos(1),pwidth,pheight],...
                             'FontSize', 8);
            plot(pfb{t}{stateOrder(s)},unit,'mean','text',mrate,'colorMap',@jet,'nanColor',nanColor);
            set(sp(end),'YTickLabels',{});
            set(sp(end),'XTickLabels',{});    
            title(regexprep(pfb{t}{stateOrder(s)}.parameters.states,'\&theta',''));
        end

        drawnow();

        % PLOT space and theta phase jointly restricted bhv rate maps ---------------------------------------
        for p = 1:(pfs{t}.adata.binSizes(2)/2),
            sp(end+1) = axes('Units','centimeters',                               ...
                             'Position',[xpos(2+p),ypos(3),pwidth,pheight],       ...
                             'FontSize', 8);
            imagescnan({pfd{t,1}.adata.bins{:},                                   ...
                        sq(mean(rmap(:,phzOrder(p)-1:phzOrder(p),:),2))'},        ...
                       [0,mrate],[],false,nanColor,'colorMap',@jet);
            axis('xy');
            axis('tight');
            set(sp(end),'YTickLabels',{});
            set(sp(end),'XTickLabels',{});        
        end 
        text(pfs{t}.adata.bins{1}(end)-0.4*diff(pfs{t}.adata.bins{1}([1,end])),   ...
             pfs{t}.adata.bins{3}(end)-0.10*diff(pfs{t}.adata.bins{3}([1,end])),  ...
             sprintf('%2.1f',mrate),                                              ...
             'Color','w','FontWeight','bold','FontSize',8)

        drawnow();

        % CREATE reference bars
        axes(fax)
        % SUPERAXIS Theta phase 
        line([sp(end-7).Position(1),sum(sp(end).Position([1,3]))],...
             sum(sp(end).Position([2,4])).*[1,1]+0.5,'LineWidth',1);
        text(mean([sp(end-7).Position(1),sum(sp(end).Position([1,3]))]),...
             sum(sp(end-7).Position([2,4]))+0.75,'Theta Phase','HorizontalAlignment','center','FontSize',8);
        % TICKS theta phase
        text(sp(end-7).Position(1),sum(sp(end-7).Position([2,4]))+0.25,...
             '0^{o}','HorizontalAlignment','center','FontSize',8);
        text(mean([sp(end-7).Position(1),sum(sp(end).Position([1,3]))]),...
             sum(sp(end-7).Position([2,4]))+0.25,'180^{o}','HorizontalAlignment','center','FontSize',8);
        text(sum(sp(end).Position([1,3])),...
             sum(sp(end-7).Position([2,4]))+0.25,'360^{o}','HorizontalAlignment','center','FontSize',8);


        cb = colorbar('Units','centimeters','Position',[sum(sp(end).Position([1,3]))+0.1,sp(end).Position(2),0.25,sp(end).Position(4)]);
        colormap(fax,'jet')
        cb.Ticks = [0,1];
        cb.TickLabels = {'min','max'};


% BSI -------------------------------------------------------------------------------------
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(5),ypos(8),pwidth*4,pheight*4],...
                         'FontSize', 8);

        imagesc(0:pfs{1}.adata.binSizes(2),-tshifts,sq(bhvInfoTS(usi,[9:16,1:8],:))');
        axis('xy')
        sp(end).XTick = [0,8,16];
        sp(end).XTickLabels = {};
        sp(end).YAxisLocation = 'right';
        ylabel('Time (ms)');
        title({'Theta Resolved and Time Shifted','Behavioral Space Information'})

        axes(fax);
        p = patch(sum(sp(end).Position([1,3]))+[-0.05,0.05,0.05,0.1,0,-0.1,-0.05]+1.5,...
                  (sp(end).Position(2)+sum(sp(end).Position([2,4])))./2+[2.5,2.5,0.35,0.35,0.15,0.35,0.35]-3,...
                  nanColor);
        uistack(p,'top')
        tx = text(sum(sp(end).Position([1,3]))+2,...
                  (sp(end).Position(2)+sum(sp(end).Position([2,4])))./2+1.5,...
                  'Future','Rotation',90,'HorizontalAlignment','center');
        uistack(tx,'top')

        p = patch(sum(sp(end).Position([1,3]))+[-0.05,0.05,0.05,0.1,0,-0.1,-0.05]+1.5,...
                  (sp(end).Position(2)+sum(sp(end).Position([2,4])))./2-[2.5,2.5,0.35,0.35,0.15,0.35,0.35]+3,...
                  nanColor);
        uistack(p,'top')
        tx = text(sum(sp(end).Position([1,3]))+2,...
                  (sp(end).Position(2)+sum(sp(end).Position([2,4])))./2-1.5,...
                  'Past','Rotation',90,'HorizontalAlignment','center');
        uistack(tx,'top')
        drawnow();



% THPP : firing rate theta phase preference -----------------------------------------------------
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(5),ypos(9),pwidth*4,pheight*1],...
                         'FontSize', 8);
        plot(pftp{t}.data.rateMap([9:16,1:8],ismember(pftp{t}.data.clu,unit)),'-+');
        sp(end).XTick = [0,8,16];
        sp(end).XTickLabels = {};
        xlim(sp(end),[0.5,16.5]);


% BAE : PLOT behavior analysis eigenvectors -----------------------------------------------------
        numComp = size(eigVec{1},2);
        pfindex = 1;

        bins = pfd{1}.adata.bins;

        % LOAD Behavioral state contours
        [stateContourMaps,stateContourHandles] =                           ...
            bhv_contours(sessionListName,                                  ... sessionListName
        'fet_HB_pitchB',                                  ... featureSet
        [1,2],                                            ... featureInd
        {linspace(-2,2,50),linspace(-2,2,50)},            ... featureBin
        'Ed05-20140529.ont.all',                          ... referenceTrial
        {'lloc+lpause&theta','hloc+hpause&theta',         ... states
         'rear&theta'},                                  ...
            'wcr'                                             ... stateColors
        );

        fpc  = cell([1,numComp]);
        for i = 1:numComp,
            fpc{i} = nan(size(validDims{pfindex}));
            fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i);
        end

        fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

        for e = 1:3,
            sp(end+1) = axes('Units','centimeters',...
                             'Position',[xpos(1),ypos(4+e),pwidth,pheight],...
                             'FontSize', 8);    
            
            imagescnan({bins{:},abs(reshape_eigen_vector(fpc{e},pfd(1,pfindex)))},...
                       fpcMinMax,'linear',false,nanColor,1,1);                % PRINT eigenvectors
            axis('xy');
            axis('tight');
            hold('on');    
            for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
                copyobj(stateContourHandles{s},sp(end));
            end
            sp(end).YTickLabel = {};
            sp(end).XTickLabel = {};
            sp(end).Color = nanColor;
            ylabel(['F',num2str(e)])
            xlim(pfd{1}.adata.bins{1}([1,end]));
            xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
            ylim(pfd{1}.adata.bins{2}([1,end]));
            if e==1, title({'Bhv Rate Map','Factor Analysis'});  end
        end
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(1),ypos(4+e+1),pwidth,pheight],...
                         'FontSize', 8);    
        plot(eigVar{1}(1:5,4),'-+');
        xlim([0,6]);
        ylabel('Variance');

% $$$ axes(fax);
% $$$ patch(mean([sp(end-2).Position(1),sum(sp(end-2).Position([1,3]))])+[-0.1,0.1,0.1,0.2,0,-0.2,-0.1],...
% $$$       sum(sp(end-2).Position([2,4]))+[1,1,0.35,0.35,0.15,0.35,0.35]-0.1,...
% $$$       nanColor);


        % BTE PLOT bhv theta resolved eigenvectors -----------------------------------------------------
        fpcMinMax = [min(nonzeros(reshape(eigVecs(usi,1:3,:,:),[],1))),...
                     max(nonzeros(reshape(eigVecs(usi,1:3,:,:),[],1)))];    
        
        for e = 1:3,
            sp(end+1) = axes('Units','centimeters',...
                             'Position',[xpos(3),ypos(4+e),pwidth,pheight],...
                             'FontSize', 8);    
            
            imagescnan({pfs{1}.adata.bins{[1,3]},sq(eigVecs(usi,e,:,:))'},...
                       fpcMinMax,'linear',false,nanColor,1,1);                % PRINT eigenvectors
            axis('xy');
            axis('tight');
% $$$     hold('on');    
% $$$     for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
% $$$         copyobj(stateContourHandles{s},sp(end));
% $$$     end
% $$$     sp(end).YTickLabel = {};
% $$$     sp(end).XTickLabel = {};
% $$$     sp(end).Color = nanColor;
% $$$     ylabel(['F',num2str(e)])
% $$$     xlim(pfd{1}.adata.bins{1}([1,end]));
% $$$     xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
% $$$     ylim(pfd{1}.adata.bins{2}([1,end]));
            if e==1, title({'Theta Resolved','Factor Analysis'});  end
        end
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(3),ypos(4+e+1),pwidth,pheight],...
                         'FontSize', 8);
        plot(eigVars(usi,:),'-+');
        xlim([0,6]);
        ylabel('Variance');




        % BTS Plot Eigenscores

        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(5),ypos(10),pwidth*4,pheight],...
                         'FontSize', 8);
        imagesc(0:pfs{1}.adata.binSizes(2),1:3,sq(eigScrs(usi,1:3,[9:16,1:8])));
        sp(end).XTick = [0,8,16];
        sp(end).XTickLabels = {'0^{o}','180^{o}','360^{o}'};
        sp(end).YTick = [1,2,3];
        sp(end).YTickLabels = {'F1','F2','F3'};
        ylabel('F Scores');
        xlabel('Theta Phase (deg)');


        FigInfo = uicontrol('Parent',hfig,...
                            'Style','text',...
                            'String',{['Unit: ',num2str(unit)],...
                                       Trial.filebase,...
                                      ['stcMode: ',Trial.stc.mode],...
                                      ['eDist:   ',num2str(Trial.nq.eDist(unit))],...
                                      ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],...
                                      ['SNR:     ',num2str(Trial.nq.SNR(unit))],...
                                      ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],...
                                      ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]...
                                     },...
                            'Units','centimeters',...
                            'Position',[xpos(12),ypos(2),6,4]);
        uistack(FigInfo,'top')        
        
        drawnow();

        FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit)];
        %print(hfig,'-depsc2',fullfile(FigDir,Trial.filebase,[FigName,'.eps']));        
        print(hfig,'-dpng',  fullfile(FigDir,Trial.filebase,[FigName,'.png']));
        
    end
end
