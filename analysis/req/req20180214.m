% req20180214 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, PlotPF, Threshold
%  Description: establish new criteria for placefield occupacy 
%               threshold
%  Bugs: NA

Trial = MTATrial.validate('er01-20110719.cof.all');
Trial = MTATrial.validate('jg05-20120317.cof.all');

dbstop at 86 in PlotPF.m

defargs = get_default_args('MjgER2016','MTAApfs','struct');
%defargs.units = [];
defargs.units = 222;
%defargs.units = 31;
%defargs.units = 61;
%defargs.units = 78;
defargs.binDims = [40,40];
defargs.SmoothingWeights = [1.25,1.25];
defargs.tag = 'test_PlotPF_new';
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs = struct2varargin(defargs);
pfn = MTAApfs(Trial,defargs{:});


% IN PlotPF.m --------------------------------------------------------------------------------
SOcc = convn(Occupancy,Smoother,'same');
SCount = convn(SpikeCount,Smoother,'same');




OccupancyNan = Occupancy;
OccupancyNan(OccupancyNan==0) = nan;
SOccNan = nanconv(Occupancy,Smoother,'edge');
SpikeCountNan = SpikeCount;
SpikeCountNan(SpikeCountNan==0) = nan;
SCountNan = nanconv(SpikeCount,Smoother,'edge');

OccThresh = 0.1.^numel(binDims);%0.06;0.12;%
                 %OccThresh = .03;%0.06;0.12;%

% Occupancy is count(xy)/xySampleRate
%    should be count(xy)/xySampleRate/binArea
% OccupancyThreshold 

20/prod(binDims)

aveSpeed = 200; %mm/s

aveSpeed/binDim(1)

figure,
gthresh = 0.05^2*5;
gtind = SOcc>(gthresh);
% SPIKECOUNT smoothed 
subplot(331);imagesc([SCount]');        set(gca,'GridColor','w');colorbar();
% OCCUPANCY smoothed
subplot(332);imagesc([SOcc]');          set(gca,'GridColor','w');colorbar();
% RATEMAP smoothed
rmap = [SCount./SOcc];
rmap(~gtind) = nan;
%subplot(333);imagescnan(rmap',[0,max(rmap(:))],[],true);  set(gca,'GridColor','w');
subplot(333);imagesc(rmap');            set(gca,'GridColor','w');colorbar();

% SPIKECOUNT raw
subplot(334);imagesc(SpikeCount');      set(gca,'GridColor','w');colorbar();
% OCCUPANCY raw
subplot(335);imagesc(Occupancy');       set(gca,'GridColor','w');colorbar();
% RATEMAP post smoothing
rmap = SpikeCount./Occupancy;           
rmap(isnan(rmap)) = 0;
rmap = convn(rmap,Smoother,'same');
rmap(~gtind) = nan;
subplot(336);imagesc(rmap');            set(gca,'GridColor','w');colorbar();

rmap = [SCountNan./SOccNan];
rmap(~gtind) = nan;
%subplot(333);imagescnan(rmap',[0,max(rmap(:))],[],true);  set(gca,'GridColor','w');
subplot(337);imagesc(rmap');            set(gca,'GridColor','w');colorbar();


%subplot(336);imagescnan(rmap',[0,max(rmap(:))],[],true);            set(gca,'GridColor','w');
subplot(338);imagesc(SOcc'>(gthresh));            set(gca,'GridColor','w');

rmap = SpikeCount./Occupancy;           
rmap(isnan(rmap)) = 0;
rmap = nanconv(rmap,Smoother,'noedge');
rmap(~gtind) = nan;
subplot(339);imagesc(rmap');            set(gca,'GridColor','w');colorbar();

ForAllSubplots(['grid(''on'');axis(''xy'')'])



% TESTING PlotPF methods


Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial.load('stc','msnn_ppsvd_raux');
units = select_placefields(Trial);

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.binDims = [10,10];
defargs.SmoothingWeights = [3,3];
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs = struct2varargin(defargs);
pf = MTAApfs(Trial,defargs{:});


defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = [92,94];
defargs.tag = 'psrmr'
defargs.binDims = [10,10];
defargs.SmoothingWeights = [3,3];
defargs.states = 'rear&theta';
defargs.overwrite = true;
defargs = struct2varargin(defargs);
pfs = MTAApfs(Trial,defargs{:});


figure,
unit = 94;
subplot(221);plot(pf ,unit,'mean',true,[],false,0.5);
subplot(222);plot(pfs,unit,'mean',true,[],false,0.5);
subplot(223);plot(pf ,unit,'mean',true,[],false,0.95);
subplot(224);plot(pfs,unit,'mean',true,[],false,0.95);


figure();
for u = units,
    plot(pf,u,'mean',true,[],false,0.99);
    waitforbuttonpress();
end

xyz = preproc_xyz(Trial,'trb');




% 20180310 



soc = Occupancy;
soc(isnan(soc)) = 0;
soc = RectFilter(soc',3,1);
soc = RectFilter(soc',3,1);
gtind = soc~=0;


ormap = SpikeCount./Occupancy;           
ormap(gtind) = 0;

rmap = SpikeCount./Occupancy;           
rmap(isnan(rmap)) = 0;
rmap = convn(rmap,Smoother,'same');

srmap = SpikeCount./Occupancy;           
srmap(isnan(srmap)) = 10.^mean(log10(srmap(:)),'omitnan');
srmap = convn(srmap,Smoother,'same');

RateMap = SCount./SOcc;
RateMap(~gtind) = nan;


rmap(~gtind)=nan;
srmap(~gtind)=nan;
RateMap(~gtind) = nan;

nx = 7;
i = 1;
cmax = [0,30];
figure;
subplot(1,nx,i);i=i+1;
imagesc(Bins{:},mean(cat(3,srmap,RateMap),3)');
%imagescnan({Bins{:},mean(cat(3,srmap,RateMap),3)'},cmax,[],true);
subplot(1,nx,i);i=i+1;
imagesc(Bins{:},srmap');

subplot(1,nx,i);i=i+1;
imagesc(Bins{:},rmap');
%imagescnan({Bins{:},rmap'},cmax,[],true);
subplot(1,nx,i);i=i+1;
imagesc(Bins{:},RateMap');
%imagescnan({Bins{:},RateMap'},cmax,[],true);
subplot(1,nx,i);i=i+1;
imagesc(SOcc');
subplot(1,nx,i);i=i+1;
imagesc(Occupancy');
subplot(1,nx,i);i=i+1;
imagesc(soc');


figure,
for u = pf.data.clu,
    clf();
    subplot(131);
    plot(pf,u,'mean',true,[],false,0.9);    
    subplot(132);
    plot(pfn,u,'mean',true,[],false,0.9);    
    subplot(133);
    imagesc(pfn.adata.bins{:},[plot(pf,u,'mean',true,[],false,0.9)-plot(pfn,u,'mean',true,[],false,0.9)]');
    colorbar();
    title(num2str(u));
    waitforbuttonpress();
end

% new PlotPF is accepted