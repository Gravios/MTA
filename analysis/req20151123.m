%Primary Goals
%Feature correction and normalization between subjects



seslist = {{'er01-20110719','hand_labeled_rev1'},...
           {'jg05-20120317','hand_labeled_rev2'},...
           {'Ed01-20140707','hand_labeled_rev1'},...
           {'Ed03-20140624','hand_labeled_rev1'}};

seslist = {...
           {'jg05-20120309','auto_wbhr'},...
           {'jg05-20120310','auto_wbhr'},...
           {'jg05-20120311','auto_wbhr'},...
           {'jg05-20120315','auto_wbhr'},...           
           {'jg05-20120317','auto_wbhr'},...
           {'jg05-20120324','auto_wbhr'},...    
           {'jg05-20120325','auto_wbhr'}};

nses = numel(seslist);

cs = jet(numel(seslist));
Trial = {}; Stc = {}; xyz = {}; fvel = {}; ang = {};
for i = 1:nses
    % Load Trials into cell arrray
    Trial{i} = MTATrial(seslist{i}{1});   
    % Load hand labeled state collections into a cell array
    Stc  {i} = Trial{i}.load('stc',seslist{i}{2});
    % Load xyz marker postions of each session into a cell array
    xyz  {i} = Trial{i}.load('xyz');    
    % Calculate lower spine marker speed from xyz position and low
    % pass filter with upper limit at 2.4Hz
    fvel {i} = xyz{i}.vel([],[1,2]).filter('ButFilter',3,2.4,'low');
    % Transform lower spine marker speed to a log scale
    fvel {i}.data(fvel{i}.data<0.01) = .01;
    fvel {i}.data = log10(fvel{i}.data);
    ang  {i} = create(MTADang,Trial{i},xyz{i});
end


% Correction of Height
mz = {};
veds = linspace(-1,2,100);
for i = 1:nses
    lind = ang{i}(:,3,4,2)<0;
    [~,ind_v_b] = histc(fvel{i}(:,'spine_lower'),veds);
    [~,ind_v_h] = histc(fvel{i}(:,'head_back'),veds);
    mind = nniz([ind_v_b,ind_v_h]);
    mz{i} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       xyz{i}(mind&lind,1,3),...
                       [numel(veds),numel(veds)],...
                       @median);
    sz{i} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       xyz{i}(mind&lind,1,3),...
                       [numel(veds),numel(veds)],...
                       @std);
end
% $$$ figure,imagesc([mz{1}]'),axis xy,colormap jet
% $$$ figure,imagesc(sz{1}'),axis xy,colormap jet
% $$$ figure,imagesc([mz{1}-mz{2}]'),axis xy,colormap jet
nnz = mz{2}~=0&mz{i}~=0;
mzd = mz{2}-mz{m};
szd = sz{2}(:)<10&sz{i}(:)<10;
% $$$ deds= linspace(-100,100,100);
% $$$ figure,bar(deds,histc(mzd(nnz(:)&szd),deds),'histc');
z_shift = median(mzd(nnz(:)&szd));


% Correction of pitch
mz = {};
veds = linspace(-1,2,100);
for i = 1:nses
    lind = ang{i}(:,3,4,2)<0;
    [~,ind_v_b] = histc(fvel{i}(:,'spine_lower'),veds);
    [~,ind_v_h] = histc(fvel{i}(:,'head_back'),veds);
    mind = nniz([ind_v_b,ind_v_h]);
    mz{i} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       ang{i}(mind&lind,5,7,2),...
                       [numel(veds),numel(veds)],...
                       @circ_median);
    sz{i} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       ang{i}(mind&lind,5,7,2),...
                       [numel(veds),numel(veds)],...
                       @circ_std);
end
% $$$ figure,imagesc([mz{2}]'),axis xy,colormap jet
% $$$ figure,imagesc(sz{2}'),axis xy,colormap jet
% $$$ figure,imagesc([mz{2}-mz{i}]'),axis xy,colormap jet
i = 1
nnz = mz{2}~=0&mz{i}~=0;
mzd = mz{2}-mz{i};
szd = sz{2}(:)<10&sz{i}(:)<.2;
deds= linspace(-pi/2,pi/2,100);
figure,bar(deds,histc(mzd(nnz(:)&szd),deds),'histc');

z_shift = median(mzd(nnz(:)&szd));


%%Start here

% Overlayed histograms of lower spine height
eds = linspace(-10,80,100);
eds = linspace(40,130,100);
eds = linspace(50,170,100);
eds = linspace(50,200,100);
figure,hold on,hout = {};
m = 1;
for i = 1:nses,
% $$$     ind = Stc{i}{'a'};
    ind = fvel{i}(:,1)>0.5&fvel{i}(:,7)>0.5;
    hout{i} = xyz{i}(ind,m,3);
    hn = bar(eds,histc(hout{i},eds)./sum(ind),'histc');%diff(ind.data,1,2)
    hn.FaceColor = cs(i,:); hn.FaceAlpha = .4; hn.EdgeAlpha = 0; % Modifyle
                                                                 % bar properties
    hout{i}(hout{i}>80|hout{i}<0) = [];
end
legend(cellfun(@subsref,Trial,repmat({substruct('.','name')},size(Trial)),'uniformoutput',false))
    

% Vertically arrange individual histograms of lower spine hight
figure,hold on
eds = linspace(-10,80,100);
for i = 1:nses,
    subplot(nses,1,i);
    ind = Stc{i}{'a'};
    bar(eds,histc(xyz{i}(ind,1,3),eds)./sum(diff(ind.data,1,2)),'histc');
    title(Trial{i}.filebase);
end



% JPDF of lower spine height vs lower spine speed over all hand
% labeled sessions 
m = 2;
out = {};
switch m
  case 1, eds = {linspace(-1,2,100),linspace(0,80,100)};
  case 2, eds = {linspace(-1,2,100),linspace(60,110,100)};
  case 3, eds = {linspace(-1,2,100),linspace(0,80,100)};
  case 4, eds = {linspace(-1,2,100),linspace(0,80,100)};
end

figure,hold on,
for i = 1:nses,
    sp = [mat2cell([floor(sqrt(numel(Trial))),ceil(sqrt(numel(Trial)))],1,ones(1,2)),i];subplot(sp{:});
    ind = Stc{i}{'w'};
    hist2([fvel{i}(ind,1),xyz{i}(ind,m,3)],eds{1},eds{2});
    xlabel('log10 speed of lower spine [log10(cm/s)]');
    ylabel('hight of lower spine [mm]');
    [out{i},xbin,ybin]= hist2([fvel{i}(ind,1),xyz{i}(ind,m,3)],eds{1},eds{2});
    title(Trial{i}.filebase);
end
suptitle('JPDF of lower spine height vs lower spine speed over all hand labeled Sessions');




% Checking how mean values of states change with state duration


figure,hold on
for s = 1:4
mz = [];
md = [];
for i = Stc{s}{'p'}.data',
    md(end+1) = log10(i(2)-i(1));
    mz(end+1) = median(xyz{s}(i(1):i(2),1,3));
end
plot(md,mz,['.',cs(s)]);
end


co = xcorr2(out{2},out{i});
figure,imagesc(co'),axis xy
mind = LocalMinimaN(-co,-1e5,1e3);
ybi = [1,mind(2)-round(size(co,1)/2)];
if ybi(2)<0,ybi = fliplr(abs(ybi));end
zshift = diff(ybin(ybi));





%Some weighted averaging of markers to estimate spine 
figure,
plot(sq(xyz{1}(1,1:4,1)),sq(xyz{1}(1,1:4,2)))
hold on
plot([mean(sq(xyz{1}(1,1:2,1)));mean(sq(xyz{1}(1,2:3,1)))],...
     [mean(sq(xyz{1}(1,1:2,2)));mean(sq(xyz{1}(1,2:3,2)))])

