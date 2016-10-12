
%% Figure 4 - Classical 2-D phase precession
Trial = MTATrial.validate('Ed10-20140817.cof.all');
Stc = Trial.load('stc','NN0317R');


Trial = 'Ed10-20140817.cof.gnd';

Trial = MTATrial.validate(Trial);
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = 'exp_phase_2d_drz';
mkdir(fullfile(OwnDir,FigDir))

% Testing Args
%Trial = MTATrial.validate('jg05-20120310');
%chans = 16;
%phase_chan = 1;
%stcMode = 'NN0317R';
%varargin = {};

[chans,phase_chan,stcMode] = DefaultArgs(varargin,{72,1,'NN0317R'});

Trial.stc = Stc;

% Load Position
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

% Load Units
pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;
% Reduce clu list based on theta pfs max rate
%units = pft.data.clu(mrt>1);

%unitt = select_units(Trial,18);
Trial.load('nq');
% $$$ units = units(Trial.nq.SNR(units)>.75);
units = find(Trial.nq.AmpSym<0&Trial.nq.SNR>.75&mrt>1)';

% Load LFP
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);
lfp.resample(xyz);
tbp_phase = lfp.phase;

% Load Stc
Trial.load('stc',stcMode);
%states = Trial.stc.list_state_attrib;
states = {'loc','rear','pause','lloc','hloc','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
nsts = numel(states);


spk = {};
pfs = {};
pmr = {};
pmp = {};
wpmr ={};
pfd = {};
DRZ = {};

% $$$ overwrite = false;
% $$$ tag = '';
% $$$ binDims = [20,20];
% $$$ smoothingWeights = [2.2,2.2];
% $$$ type = 'xy';

for s = 1:nsts,        
    % Load/Create place fields

    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{s};
    defargs = struct2varargin(defargs);        
    pfs{s} = MTAAknnpfs_bs(Trial,defargs{:});      


% $$$     pfs{s} = MTAApfs(Trial,...
% $$$                      units,...
% $$$                      states{s},...
% $$$                      overwrite,...
% $$$                      tag,...
% $$$                      binDims,...
% $$$                      smoothingWeights,...
% $$$                      type);
    
    % Load spikes 
    spk{s} = Trial.spk.copy;
    spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');

    
    % Get the mean firing rate for each xy position along trajectory 
    wpmr{s} = zeros(xyz.size(1),numel(units));

    [~,indx] = min(abs( repmat(pfs{s}.adata.bins{1}',xyz.size(1),1)...
                        -repmat(xyz(:,Trial.trackingMarker,1),1,numel(pfs{s}.adata.bins{1}))),...
                   [],2);

    [~,indy] = min(abs( repmat(pfs{s}.adata.bins{2}',xyz.size(1),1)...
                        -repmat(xyz(:,Trial.trackingMarker,2),1,numel(pfs{s}.adata.bins{2}))),...
                   [],2);

    rateMapIndex = sub2ind(pfs{s}.adata.binSizes',indx,indy);
    for unit = units,
        rateMap = pfs{s}.plot(unit,'mean');      %  for MTAApfs
                                          %rateMap = rot90(rot90(rateMap)'); % for MTAAknnpfs
        wpmr{s}(:,unit==units) = rateMap(rateMapIndex);
    end

    % Get the peak firing rate and position of each place field
    [pmr{s},pmp{s}] = pfs{s}.maxRate(units,'mean');
    pmr{s} = repmat(pmr{s}(:)',xyz.size(1),1);


    % Get the rat's heading 
    pfds = [];
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
    end

    pfd{s} = zeros(size(pfds));
    pfd{s}(abs(pfds)<=pi/2)=-1;
    pfd{s}(abs(pfds)>pi/2)=1;

    % Calculate DRZ 
    DRZ{s} = pfd{s}.*(1-wpmr{s}./pmr{s});

end



width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;


% 2d DRZ figure
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

aIncr = false;
hfig = figure(38384);
hfig.Units = 'centimeters';
hfig.Position = [1,1,40,24];    
hfig.PaperPositionMode = 'auto';

unit = units(1);
while unit~=-1,
    
    clf
    for s = 1:nsts,
        res = spk{s}(unit);
        
        if numel(res) <50,continue,end
        res(res>xyz.size(1))=[];            
        drzspk = DRZ{s}(res,unit==units);
        phzspk = tbp_phase(res,phase_chan);
        
        gind = ~isnan(drzspk)&~isnan(phzspk);
        
        subplot2(9,nsts,[1,2],s);
        plot(xyz(res,Trial.trackingMarker,1),xyz(res,Trial.trackingMarker,2),'.b');
        xlim([-500,500]),ylim([-500,500])
        title(states{s})
        
        subplot2(9,nsts,[4,5],s);
        pfs{s}.plot(unit,'mean','isCircular',true);
        text(pfs{s}.adata.bins{1}(end)-250,pfs{s}.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)

        hold on,plot(pmp{s}(unit==units,1),pmp{s}(unit==units,2),'w*')
        title(num2str(unit))
        
        if sum(gind)>10,
            subplot2(9,nsts,[7,9],s);
            hold('on');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind))+720,'b.');
            xlim([-1,1]),
            ylim([-180,900])
        end
    end
    
    
    FigName = [Trial.filebase,'_',FigDir,'_unit',num2str(unit)];
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));

    unit = figure_controls(hfig,unit,units,aIncr);        
    
end

