
%% Figure 4 - Classical 2-D phase precession

Trial = MTATrial.validate('jg05-20120309.cof.all');
chans = 63;

Trial = MTATrial.validate('jg05-20120310.cof.all');
chans = 63;

Trial = MTATrial.validate('jg05-20120311.cof.all');
chans = 63;

Trial = MTATrial.validate('jg05-20120317.cof.all');
chans = 50;

Trial = MTATrial.validate('ER06-20130613.cof.all');
chans = 32;

Trial = MTATrial.validate('ER06-20130614.cof.gnd');
chans = 32;

Trial = MTATrial.validate('Ed10-20140817.cof.gnd');
chans = [8,16,17,25];

%Trial = MTATrial.validate(Trial);
stcMode = 'msnn_ppsvd_raux';
Trial.load('stc',stcMode);
Stc = Trial.stc.copy();
OwnDir = '/storage/gravio/nextcloud/MjgER2016/figures/figure3';
FigDir = 'exp_phase_2d_drz_theta';
mkdir(fullfile(OwnDir,FigDir))
phase_chan = 1;

% Testing Args
%Trial = MTATrial.validate('jg05-20120310');
%chans = 68;
%varargin = {};
%[chans,phase_chan,stcMode] = DefaultArgs(varargin,{[66,72,78,84],1,'NN0317R'});

% Load Stc
states = {'rear&theta','loc&theta','lloc&theta','hloc&theta',...
          'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};
nsts = numel(states);

% Load Position
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

% Load Units
pft = pfs_2d_theta(Trial);
[mrt,pmp] = pft.maxRate;

pfstats = compute_pfstats_bs(Trial);
patchCOM = [];
for unit = units,
    patchCOM(end+1,:) = fliplr(sq(mean(pfstats.peakPatchCOM(8,:,pfstats.cluMap==unit,:)))');
end

% Reduce clu list based on theta pfs max rate
units = select_units(Trial,18);
units = units(mrt(units)>1);

% Load LFP
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);
lfp.resample(xyz);
tbp_phase = lfp.phase([6,12]);


for s = 1:nsts,        
    % Load spikes 
    spk{s} = Trial.spk.copy;
    spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');
end

% COMPUTE directed rate zones
drz = compute_drz(Trial,pft,units);
ddz = compute_ddz(Trial,pft,units);


% load autocorrelogram
[accg,tbins] = autoccg(MTASession.validate(Trial.filebase));


% FIGURE plot place field, 2d DRZ vs theta phase
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

aIncr = true;
hfig = figure(38384);
hfig.Units = 'centimeters';
hfig.Position = [1,1,40,20];    
hfig.PaperPositionMode = 'auto';

unit = units(1);
distThresh = 250;
while unit~=-1,
 
   clf
    for s = 1:nsts,
        res = spk{s}(unit);
        % Plot unit auto correlogram
        if s == 1,
            subplot2(9,nsts,[9],s);
            bar(tbins,accg(:,unit));axis tight;            
        end

        if numel(res) <50,continue,end
        res(res>xyz.size(1))=[];

        drzspk = drz(res,unit==units);
        ddzspk = ddz(res,unit==units);        
        phzspk = tbp_phase(res,phase_chan);
        gind = ~isnan(drzspk)&~isnan(phzspk);

        subplot2(9,nsts,[1,2],s);  hold('on');
        plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.m');
        xlim([-500,500]),ylim([-500,500])
        title(states{s})

        subplot2(9,nsts,[4,5],s);  hold('on');
        pft.plot(unit,'mean',[],mrt(unit==pft.data.clu).*1.5,'isCircular',true);
        %plot(pmp(unit==pft.data.clu,1),pmp(unit==pft.data.clu,2),'w*')
        plot(patchCOM(unit==units,1),patchCOM(unit==units,2),'w*')
        title(num2str(unit))

        % Plot phase drz relationship
        if sum(gind)>10,
            subplot2(9,nsts,[7,8],s);
            hold('on');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
            xlim([-1,1]),
            ylim([-180,540])
        end
% $$$         if sum(gind)>10,
% $$$             subplot2(9,nsts,[7,8],s);
% $$$             hold('on');
% $$$             plot(ddzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
% $$$             plot(ddzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
% $$$             xlim([-600,600]),
% $$$             ylim([-180,540])
% $$$         end

    end
    
    FigName = [Trial.filebase,'_',FigDir,'_nearPatchCenter',num2str(distThresh),'_unit',num2str(unit)];
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    unit = figure_controls(hfig,unit,units,aIncr);        
end




% FIGURE plot place field, 2d DRZ vs theta phase
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

aIncr = true;
hfig = figure(38384);
hfig.Units = 'centimeters';
hfig.Position = [1,1,40,20];    
hfig.PaperPositionMode = 'auto';

unit = units(1);
distThresh = 250;
while unit~=-1,
 
   clf
    for s = 1:nsts,
        res = spk{s}(unit);
        % Plot unit auto correlogram
        if s == 1,
            subplot2(9,nsts,[9],s);
            bar(tbins,accg(:,unit));axis tight;            
        end

        if numel(res) <50,continue,end
        res(res>xyz.size(1))=[];
% $$$         ares=res;
% $$$         ares(pfdist(ares,unit==units)<distThresh)=[];
% $$$         res(pfdist(res,unit==units)>=distThresh)=[];            

        drzspk = drz(res,unit==units);
        phzspk = tbp_phase(res,phase_chan);
        gind = ~isnan(drzspk)&~isnan(phzspk);

        subplot2(9,nsts,[1,2],s);hold on
% $$$         plot(xyz(ares,Trial.trackingMarker,1),xyz(ares,Trial.trackingMarker,2),'.b');
        plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.m');
        xlim([-500,500]),ylim([-500,500])
        title(states{s})

        subplot2(9,nsts,[4,5],s);
        pft.plot(unit,'mean',[],mrt(unit==pft.data.clu).*1.5,'isCircular',true);
        hold on,plot(pmp(unit==pft.data.clu,1),pmp(unit==pft.data.clu,2),'w*')
        title(num2str(unit))

        % Plot phase drz relationship
        if sum(gind)>10,
            subplot2(9,nsts,[7,8],s);
            hold('on');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
            xlim([-1,1]),
            ylim([-180,540])
        end
    end
    
    FigName = [Trial.filebase,'_',FigDir,'_nearPatchCenter',num2str(distThresh),'_unit',num2str(unit)];
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    unit = figure_controls(hfig,unit,units,aIncr);        
end





%Phase prec
pfs = {};
for s = 1:nsts,        
    % Load/Create place fields
    defargs = get_default_args('MjgEdER2016','MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{s};
    defargs.numIter = 1001;
    defargs = struct2varargin(defargs);        
    pfs{s} = MTAAknnpfs_bs(Trial,defargs{:});      
end

% use unit 33 from Ed10-20140817.cof.gnd
distThresh = 250;
mResults   = 3;
fitRange   = -pi:0.01:pi;
P          = nan([numel(units),nsts,mResults,2]);
phzStats   = nan([numel(units),nsts,mResults,2]);



if reportFig, 
    hfig = figure(20161028);
    hfig.Units    = 'centimeters';
    hfig.Position = [0,0,11.5,25];
    hfig.PaperPositionMode = 'auto';
    autoincr = true;

end

unit = units(1);
while unit~=-1,
    if reportFig,clf(hfig),end
    for s = 1:nsts

        res = spk{s}(unit);
        res(res>xyz.size(1))=[];
        ares=res;
        ares(pfdist(ares,unit==units)< distThresh) = [];
        res (pfdist(res, unit==units)>=distThresh) = [];            

        drzspk = DRZ(res,unit==units);
        phzspk = tbp_phase(res,phase_chan);
        gind = ~isnan(drzspk)&~isnan(phzspk);

        if numel(gind)<20,
            continue;
        end
        
        lin = drzspk(gind);
        circ = phzspk(gind);
        x=[min(lin),max(lin)];

        cosPart = sum(cos(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
        sinPart = sum(sin(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
        R = sqrt((cosPart./length(circ)).^2+...
                 (sinPart./length(circ)).^2 );


        [lmi,lmv] = LocalMinima(-R',0,0,mResults);
        lmid = find(~isnan(lmi));
        % P:  Regression Parm: [Slope         ,Offset                          ]        
        P(unit==units,s,lmid,:) = [fitRange(lmi(lmid))',atan2(sinPart(lmi(lmid)),cosPart(lmi(lmid)))'];

        % Collect residuals of the theta model for each state
        phi = nan([numel(lin),mResults]);
        for r = lmid',
            [phi(:,r)] = polyval([2*pi;1].*sq(P(unit==units,1,r,:)),lin);
        end
        residuals = bsxfun(@minus,circ,phi);
        residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
        residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;

        if reportFig, 
            subplot2(nsts,4,s,1); pft.plot(unit,'mean',[],mrt(unit).*1.5,'isCircular',true);
                                  ylabel(states{s});
            subplot2(nsts,4,s,2); 
                      plot(fitRange,R);   xlim([-pi,pi]); hold on
                      LocalMinima(-R',numel(fitRange),0,0,mResults);
                      for r = lmid',plot(fitRange(lmi(r)),-lmv(r),'o'),end
            subplot2(nsts,4,s,3); plot(lin,circ,'.'); ylim([-pi,pi]);
                      hold('on'); 
                      for r = lmid',plot(x,2*pi*P(unit==units,s,r,1)*x+P(unit==units,s,r,2),'-'),end
            subplot2(nsts,4,s,4); hist(residuals(:,1),100);xlim([-pi,pi]);

        end

        phzStats(unit==units,s,:,:) = [circ_mean(residuals);circ_std(residuals)]';
    end
    suptitle([Trial.filebase,' unit: ',num2str(unit)]);

    FigName = [Trial.filebase,'-phase_prcs-unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end




