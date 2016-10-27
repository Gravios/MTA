
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
chans = 16;

%Trial = MTATrial.validate(Trial);
stcMode = 'NN0317R';
Stc = Trial.load('stc',stcMode);
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = 'exp_phase_2d_drz_theta';
mkdir(fullfile(OwnDir,FigDir))
phase_chan = 1;

% Testing Args
%Trial = MTATrial.validate('jg05-20120310');
%chans = 68;
%varargin = {};
%[chans,phase_chan,stcMode] = DefaultArgs(varargin,{[66,72,78,84],1,'NN0317R'});

% Load Stc
states = {'rear','loc','lloc','hloc','pause','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states = cat(2,{'theta-sit-groom'},states,{'pause-theta'});
nsts = numel(states);

% Load Position
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

% Load Units
pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;

% Reduce clu list based on theta pfs max rate
units = select_units(Trial,18);
units = units(mrt(units)>1);

% Load LFP
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);
lfp.resample(xyz);
tbp_phase = lfp.phase([6,12]);



spk = {};
pfdist = [];
pmr  = [];
pmp  = [];
wpmr = [];
pfd  = [];
DRZ  = [];


for s = 1:nsts,        
    % Load spikes 
    spk{s} = Trial.spk.copy;
    spk{s}.create(Trial,xyz.sampleRate,states{s},[],'deburst');
end


    
% Get the mean firing rate for each xy position along trajectory 
wpmr = zeros(xyz.size(1),numel(units));

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
    wpmr(:,unit==units) = rateMap(rateMapIndex);
end

    % Get the peak firing rate and position of each place field
[pmr,pmp] = pft.maxRate(units,'mean');
pmr = repmat(pmr(:)',xyz.size(1),1);


% Get the rat's 
pfds = [];
pfdd = [];
for unit = units
    pfhxy = xyz(:,{'head_back','head_front'},:);
    pfhxy = cat(2,pfhxy,permute(repmat([pmp(unit==units,:),0],xyz.size(1),1),[1,3,2]));
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

pfdist = pfdd;
pfd = zeros(size(pfds));
pfd(abs(pfds)<=pi/2)=-1;
pfd(abs(pfds)>pi/2)=1;

% Calculate DRZ 
DRZ = pfd.*(1-wpmr./pmr);



%spk,DRZ,pfdist

% load autocorrelogram
[accg,tbins] = autoccg(MTASession.validate(Trial.filebase));


% FIGURE 2d DRZ 
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

aIncr = true;
hfig = figure(38384);
hfig.Units = 'centimeters';
hfig.Position = [1,1,55,24];    
hfig.PaperPositionMode = 'auto';

unit = units(1);
distThresh = 350;
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
        ares(pfdist(ares,unit==units)<distThresh)=[];
        res(pfdist(res,unit==units)>=distThresh)=[];            

        drzspk = DRZ(res,unit==units);
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
        hold on,plot(pmp(unit==units,1),pmp(unit==units,2),'w*')
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
pause(0.2)    
end




mrm = pft.plot(unit,'mean',[],mrt(unit).*1.5,'isCircular',true);
srm = pft.plot(unit, 'std',[],mrt(unit).*1.5,'isCircular',true);



% use unit 33 from Ed10-20140817.cof.gnd

lin = drzspk(gind);
circ = phzspk(gind);

% $$$ a=-pi:0.001:pi;
% $$$ cosinepart=zeros(length(a),1);
% $$$ sinepart=zeros(length(a),1);
% $$$ R=zeros(length(a),1);
% $$$ for i=1:length(a)
% $$$     cosinepart(i)=sum(cos(circ-(2*pi*a(i)*lin)));
% $$$     sinepart(i)=sum(sin(circ-(2*pi*a(i)*lin)));
% $$$     firstterm=(cosinepart(i)/length(circ))^2;
% $$$     secondterm=(sinepart(i)/length(circ))^2;
% $$$     R(i)=sqrt(firstterm+secondterm);
% $$$ end
% $$$ figure
% $$$ plot(a,R);
% $$$ slope=a(R==max(R));
% $$$ offset=atan2(sinepart(R==max(R)),cosinepart(R==max(R)));
% $$$ x=min(lin):0.1:max(lin);


figure,hold on
plot(lin,circ,'.')
plot(x,2*pi*slope*x+offset,'r-')


[P] = polyval([2*pi*slope,offset],lin);
residuals = circ-P;
residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;

figure,hist(residuals,100)



