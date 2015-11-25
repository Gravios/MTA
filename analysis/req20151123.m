%Primary Goals
%Feature correction and normalization between subjects



seslist = {{'er01-20110719','hand_labeled_rev1'},...
           {'jg05-20120317','hand_labeled_rev2'},...
           {'Ed01-20140707','hand_labeled_rev1'},...
           {'Ed03-20140624','hand_labeled_rev1'}};
nses = numel(seslist);

cs = 'rbgc';
Trial = {}; Stc = {}; xyz = {};
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
    fvel{i}(fvel{i}<0.01) = .01;
    fvel {i}.data = log10(fvel{i}.data);
end


% Overlayed histograms of lower spine height
eds = linspace(0,80,100);
figure,hold on
for i = 1:nses,
    ind = Stc{i}{'a'};
    hn = bar(eds,histc(xyz{i}(ind,1,3),eds)./sum(diff(ind.data,1,2)),'histc');
    hn.FaceColor = cs(i); hn.FaceAlpha = .4; hn.EdgeAlpha = 0; % Modifyle
                                                               % bar properties
end
legend(cellfun(@subsref,Trial,repmat({substruct('.','name')},size(Trial)),'uniformoutput',false))
    

% Vertically arrange individual histograms of lower spine hight
figure,hold on
eds = linspace(0,80,100);
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

