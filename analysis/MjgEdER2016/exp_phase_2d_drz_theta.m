
%% Figure 4 - Classical 2-D phase precession
Trial = MTATrial.validate('Ed10-20140817.cof.gnd');
Stc = Trial.load('stc','NN0317R');

Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial = MTATrial.validate('jg05-20120309.cof.all');

Trial = MTATrial.validate(Trial);
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = 'exp_phase_2d_drz_theta';
mkdir(fullfile(OwnDir,FigDir))

% Testing Args
%Trial = MTATrial.validate('jg05-20120310');
%chans = 68;
%phase_chan = 1;
%stcMode = 'NN0317R';
%varargin = {};

[chans,phase_chan,stcMode] = DefaultArgs(varargin,{72,1,'NN0317R'});

% Load Stc
Trial.load('stc',stcMode);
%states = Trial.stc.list_state_attrib;
states = {'rear','loc','lloc','hloc','pause','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states = cat(2,{'theta-sit-groom'},states,{'pause-theta'});
nsts = numel(states);

% Load Position
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

% Load Units
pft = pfs_2d_theta(Trial,[],[],true);
mrt = pft.maxRate;

% Reduce clu list based on theta pfs max rate
units = select_units(Trial,18);
units = units(mrt(units)>1);

% Load LFP
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);
lfp.resample(xyz);
tbp_phase = lfp.phase;



spk = {};
pfdist = {};
pmr = {};
pmp = {};
wpmr ={};
pfd = {};
DRZ = {};


for s = 1:nsts,        
    % Load/Create place fields

    % Load spikes 
    spk{s} = Trial.spk.copy;
    spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');

    
    % Get the mean firing rate for each xy position along trajectory 
    wpmr{s} = zeros(xyz.size(1),numel(units));

    [~,indx] = min(abs( repmat(pft.adata.bins{1}',xyz.size(1),1)...
                        -repmat(xyz(:,Trial.trackingMarker,1),1,numel(pft.adata.bins{1}))),...
                   [],2);

    [~,indy] = min(abs( repmat(pft.adata.bins{2}',xyz.size(1),1)...
                        -repmat(xyz(:,Trial.trackingMarker,2),1,numel(pft.adata.bins{2}))),...
                   [],2);

    rateMapIndex = sub2ind(pft.adata.binSizes',indx,indy);
    for unit = units,
        rateMap = pft.plot(unit,'mean');      %  for MTAApfs
                                          %rateMap = rot90(rot90(rateMap)'); % for MTAAknnpfs
        wpmr{s}(:,unit==units) = rateMap(rateMapIndex);
    end

    % Get the peak firing rate and position of each place field
    [pmr{s},pmp{s}] = pft.maxRate(units,'mean');
    pmr{s} = repmat(pmr{s}(:)',xyz.size(1),1);


    % Get the rat's heading 
    pfds = [];
    pfdd = [];
    for unit = units
        pfhxy = xyz(:,{'head_back','head_front'},:);
        pfhxy = cat(2,pfhxy,permute(repmat([pmp{s}(unit==units,:),0],xyz.size(1),1),[1,3,2]));
        pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
        
        cor = cell(1,3);
        [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
        cor = cell2mat(cor);
        
        por = cell(1,3);
        [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
        por = cell2mat(por);
        
        pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
        pfdd(:,unit==units) = por(:,3);
        
    end

    pfdist{s} = pfdd;
    pfd{s} = zeros(size(pfds));
    pfd{s}(abs(pfds)<=pi/2)=-1;
    pfd{s}(abs(pfds)>pi/2)=1;

    % Calculate DRZ 
    DRZ{s} = pfd{s}.*(1-wpmr{s}./pmr{s});

end

[accg,tbins] = autoccg(MTASession.validate(Trial.filebase));

width = pft.adata.binSizes(1);
height = pft.adata.binSizes(2);
radius = round(pft.adata.binSizes(1)/2)-find(pft.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;


% 2d DRZ figure
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

aIncr = true;
hfig = figure(38384);
hfig.Units = 'centimeters';
hfig.Position = [1,1,55,24];    
hfig.PaperPositionMode = 'auto';

unit = units(1);
distThresh = 200;
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
        ares=res;
        ares(pfdist{s}(ares,unit==units)<distThresh)=[];
        res(pfdist{s}(res,unit==units)>=distThresh)=[];            

        drzspk = DRZ{s}(res,unit==units);
        phzspk = tbp_phase(res,phase_chan);
        
        gind = ~isnan(drzspk)&~isnan(phzspk);
        
        subplot2(9,nsts,[1,2],s);hold on
        plot(xyz(ares,Trial.trackingMarker,1),xyz(ares,Trial.trackingMarker,2),'.b');
        plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.m');
        xlim([-500,500]),ylim([-500,500])
        title(states{s})
        
        subplot2(9,nsts,[4,5],s);
        pft.plot(unit,'mean',[],mrt(unit).*1.5,'isCircular',true);
        text(pft.adata.bins{1}(end)-250,pft.adata.bins{2}(end)-50,...
             sprintf('%2.1f',pft.maxRate(unit)),'Color','w','FontWeight','bold','FontSize',10)
        hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
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




mrm = pft.plot(unit,'mean',[],mrt(unit).*1.5,'isCircular',true);
srm = pft.plot(unit, 'std',[],mrt(unit).*1.5,'isCircular',true);