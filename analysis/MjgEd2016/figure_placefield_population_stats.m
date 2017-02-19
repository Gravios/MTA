
% LOAD DATA ---------------------------------------------------------------------------
% Summary
%   N := number of units
%   S := number of states
%
%   map    - Nx4 matrix with unit identification
%   nq     - Nx1 struct with array fields regarding unit characteristics   
%   mpfs   - SxN struct with array fields regarding place field characteristics
%   states - 1xS cell array of strings with behavior state labels
%            corresponding to all other varialbles with an S dim

load_pfs_unit_stats;

%--------------------------------------------------------------------------------------




% FIGURES -----------------------------------------------------------------------------
% Plotting some Group stats
% 
% Rear x Loc
% 
% Rear x Pause
%
% Loc  x Pause
%

% mean patch firing rate 


uind = nq.eDist>20&nq.Refrac<5e-4;

%uind = ':';
figure, hold on
scatter(mpfs.patchMFR(s1,uind,1),... %.*log10(mpfs.patchArea(s1,uind,1)),...
        mpfs.patchMFR(s2,uind,1),... %.*log10(mpfs.patchArea(s2,uind,1)),...
        nq.eDist(uind)/2,...
        cat(2,MakeUniformDistr(nq.Refrac(uind))/max(nq.Refrac(uind)).*double(nq.Refrac(uind)~=0),zeros([sum(uind),2])),...
        'filled');

plot([0,100],[0,100],'-m')
plot([0,50],[0,100],'-r')
plot([0,100],[0,50],'-r')
plot([0,25],[0,100],'-r')
plot([0,100],[0,25],'-r')
daspect([1,1,1])
xlim([0,25])
ylim([0,25])


s1 = 1;
s2 = 2;
uind = nq.eDist>25;
uind = nq.Refrac<5e-4;
figure, hold on
plot(mpfs.patchMFR(s1,uind,1),...
     mpfs.patchMFR(s2,uind,1),...
     '.')
% $$$ plot(mpfs.patchPFR(s1,uind,1),...
% $$$      mpfs.patchPFR(s2,uind,1),...
% $$$      '.')
plot([0,40],[0,40],'-m')
plot([0,20],[0,40],'-r')
plot([0,40],[0,20],'-r')
plot([0,10],[0,40],'-r')
plot([0,40],[0,10],'-r')
daspect([1,1,1])


s1 = 1;
s2 = 2;
figure,plot(mpfs.patchMFR(s1,:,1),...
            mpfs.patchMFR(s2,:,1),...
            '.')


s1 = 1;
s2 = 2;
figure,plot(log10(mpfs.peakFR(s1,uind,1)),...
            log10(mpfs.peakFR(s2,uind,1)),...
            '.')
daspect([1,1,1])
hold on,plot(log10(linspace(0.1,40,100)),log10(linspace(0.1,40,100)),'-m')
hold on,plot(log10(linspace(0.1,20,100)),log10(linspace(0.1,40,100)),'-r')
hold on,plot(log10(linspace(0.1,40,100)),log10(linspace(0.1,20,100)),'-r')


figure,hist(log10(mpfs.peakFR(s1,uind,1)./mpfs.peakFR(s2,uind,1)),100)


% log 10 mean patch rate
%uind = true([numel(nq.eDist),1]);
%uind = nq.eDist>30&nq.Refrac<5e-4;
uind = nq.Refrac<4e-4;
s1 = 1;
s2 = 2;
figure,
scatter(log10(mpfs.patchPFR(s1,uind,1)),... %.*log10(mpfs.patchArea(s1,uind,1)),...
        log10(mpfs.patchPFR(s2,uind,1)),... %.*log10(mpfs.patchArea(s2,uind,1)),...
        nq.eDist(uind)/2,...
        cat(2,MakeUniformDistr(nq.Refrac(uind))/max(nq.Refrac(uind)).*double(nq.Refrac(uind)~=0),zeros([sum(uind),2])),...
        'filled');
daspect([1,1,1])
hold on,plot(log10(linspace(0.1,40,1000)),log10(linspace(0.1,40,1000)),'-b')
hold on,plot(log10(linspace(0.1,20,1000)),log10(linspace(0.1,40,1000)),'-r')
hold on,plot(log10(linspace(0.1,40,1000)),log10(linspace(0.1,20,1000)),'-r')
hold on,plot(log10(linspace(0.1,20,1000)),log10(linspace(0.1,80,1000)),'-m')
hold on,plot(log10(linspace(0.1,80,1000)),log10(linspace(0.1,20,1000)),'-m')
xlim([-1.5,1.75])
ylim([-1.5,1.75])
daspect([1,1,1])
xlabel(['mean patch firing rate (log10(Hz)) during ',states{s1}, ' behavior'])
ylabel(['mean patch firing rate (log10(Hz)) during ',states{s2}, ' behavior'])



figure,hist(log10(mpfs.peakFR(s1,uind,1)./mpfs.peakFR(s2,uind,1)),100)


%--------------------------------------------------------------------------------------
