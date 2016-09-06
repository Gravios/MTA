

Slist = SessionList('ER06');

% $$$ QuickSessionSetup(Slist(4));
% $$$ 
s = MTASession.validate(Slist(4));
% $$$ pXY(s);
% $$$ pZ(s);
% $$$ 
% $$$ NeuronQuality(s);

xyz = s.load('xyz');
% $$$ xyz.data(:,:,1) = xyz.data(:,:,1)+Slist(2).xOffSet;
% $$$ xyz.data(:,:,2) = xyz.data(:,:,2)+Slist(2).yOffSet;
% $$$ xyz.save;

Tlist = SessionList('exp_flight');

QuickTrialSetup(Tlist,'overwrite',true);

Trial = MTATrial.validate('ER06-20130614.cof.gnd');
Trial.maze.boundaries(end) = 450;

display = true;
overwrite = true;
units = select_units(Trial,18,'pyr');


Trial = labelBhv_NN(Trial,'NN0317');
Trial = labelTheta(Trial,[],73);





pfs = {};
pfs{1} = MTAApfs(Trial,units,'theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);


Trial = MTATrial.validate('ER06-20130614.cof.fly');
Trial.maze.boundaries(end) = 450;

if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end
if isempty(Trial.stc.gsi('f')),Trial = labelFlight(Trial);end
%Trial = labelFlight(Trial,[],1,'set');



pfs{2} = MTAApfs(Trial,units,'theta-flight',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);

pfs{3} = MTAApfs(Trial,units,'flight&theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);







%% Plot placefields across height
slices = 1:2:19;

spOpts.width  = 2;
spOpts.height = 2;
spOpts.ny = 3;
spOpts.nx = numel(slices);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 0;
figOpts.position = [1,1,(spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding,...
                     (spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2)];

width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;
    
    

unit = units(128);
hfig = figure(393929);clf
for i  = 1:3
pf = pfs{i};
ratemap = pf.plot(unit,'isCircular',false);
ratemap(isnan(ratemap)) = -1;
for s = 1:numel(slices)
    sp(i,s) = axes('Units',spOpts.units,...
                   'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2),...
                        (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                        spOpts.width,...
                        spOpts.height]...
                   );
    hold('on')
    imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap(:,:,slices(s)).*mask');    
    axis xy
    colormap([0,0,0;parula]);
    caxis([-1,max(ratemap(:).*reshape(repmat(mask,[1,1,size(ratemap,3)]),[],1))]);
    title(num2str(round(pf.adata.bins{3}(slices(s)))))
    text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
         sprintf('%2.1f',max(max(ratemap(:,:,slices(s))))),'Color','w','FontWeight','bold','FontSize',10)
end
end





if display,    
    set(0,'defaultAxesFontSize',8,...
          'defaultTextFontSize',8)

    autoincr = false;
    hfig = figure(38384);
    hfig.Units = 'centimeters';
    hfig.Position = [1,1,40,24];    
    hfig.PaperPositionMode = 'auto';

    while unit~=-1,
        sp = [];
        for i=1:numel(pfs),
            mrxz(i) =  max(max(max(pfs{i}.plot(unit,'xz',1,'isCircular',false))));
        end

        Nmrxz = max(mrxz);
        for i=1:numel(pfs),
            subplot2(numel(pfs),2,i,1);
            pfs{i}.plot(unit,'xz',1,[0,mrxz],'isCircular',false); 
            title([pfs{i}.session.trialName,':',pfs{i}.parameters.states,': ',num2str(unit)]);
            subplot2(numel(pfs),2,i,2);
            pfs{i}.plot(unit,'xy',1,[0,mrxz],'isCircular',false); 
        end
        saveas(gcf,fullfile('/gpfs01/sirota/home/gravio/',...
                            'figures','SFN2014',...
                            ['pfsFly2_' Trial.filebase '-' num2str(unit) '.png']),'png');

        unit = figure_controls(hfig,unit,units,autoincr);
    end

end
