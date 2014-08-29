
%% Section 1 Plot placefield and unit accgs to survey the unit population
Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120317');
MTAstartup('cin','cin');Trial = MTATrial('Ed10-20140812');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;

linkSession(Trial.name,'/gpfs01/sirota/homes/eduardo/data/xyz','/gpfs01/sirota/homes/eduardo/data/rawnlx')

if isempty(Trial.stc{'t'}),
    Trial = labelAuxBhv(Trial,Trial.stc,'overwrite',true);
end


units = select_units(Trial,18,'pyr');

% $$$ Trial.stc.states{Trial.stc.gsi('c')}.data = Trial.stc{'c'}(diff(Trial.stc{'c'}.data(:,:),1,2)>20,:);
% $$$ Trial.stc.states{Trial.stc.gsi('p')}.data = Trial.stc{'p'}(diff(Trial.stc{'p'}.data(:,:),1,2)>20,:);

stss = 'twcpr';
pfs = {};
ow =true;
for s = 1:numel(stss),
    pfs{s} = MTAApfs(Trial,units,Trial.stc{stss(s)},ow);
end

[accg,tbin] = autoccg(Trial);


nx = 4;ny = 2;        
hfig = figure(3747473);
i = units(1);
while i ~=-1,
    clf
    p = 1;
    for s=1:nx*ny
        subplot(ny,nx,s);
        if s == 5,
            bar(tbin,accg(:,i)),axis tight,title(num2str(i))
        elseif p<=numel(pfs),
            pfs{p}.plot(i,[],1);
            title(pfs{p}.parameters.states)
            p = p+1;
        end
    end
    i = figure_controls(hfig,i,units);
end


%% Section 1.2 Load old low walk and high walk periods as defined
%% by hmm segmentation of head height, angle and velocity
ds = load('/gpfs01/sirota/data/bachdata/data/gravio/stasis/jg05-20120317/jg05-20120317.cof.all.stc.auto_wbhr.mat');
stsind=6;
old_stc_origin = ds.Stc.states{stsind}.origin;
hwper = ds.Stc.states{stsind}.copy;
hwper.origin = Trial.stc.states{1}.origin;
hwper.path = Trial.stc.states{1}.path;
hwper.sync = Trial.stc.states{1}.sync.copy;
thwdata = hwper.data - (Trial.stc.states{1}.origin*Trial.xyz.sampleRate-old_stc_origin);
szero = find(thwdata(:,2)<0,1,'last');
if ~isempty(szero),
    thwdata(1:szero,:) = [];
end
if thwdata(1)<=0,thwdata(1)=1;end
hwper.data = thwdata;

Trial.stc.states{end+1} = hwper;


%% Section 1.3 Various Distributions compared between states
xyz = Trial.xyz.copy;xyz.load(Trial);
vel = xyz.vel(7,[1,2])
vel.data = clip(vel.data,0,100);
figure,plot(vel.data),
Lines(Trial.stc{'p'}(:)+1,[],'c');
Lines(Trial.stc{'w'}(:),[],'g');
Lines(Trial.stc{'l'}(:)+1,[],'k');

Lines(Trial.stc{'g'}(:)+1,[],'r');
Lines(Trial.stc{'c'}(:)+1,[],'m');


figure,hist(ang(Trial.stc{'w'},5,7,2),linspace(-1.4,1.4,500));

figure,
subplot(211),hist(tan(ang(Trial.stc{'l'},5,7,2)),linspace(-1.4,1.4,500));
subplot(212),hist(tan(ang(Trial.stc{'g'},5,7,2)),linspace(-1.4,1.4,500));

figure,
subplot(211),hist(tan(ang(Trial.stc{'p'},5,7,2)),linspace(-1.4,1.4,500));
subplot(212),hist(tan(ang(Trial.stc{'c'},5,7,2)),linspace(-1.4,1.4,500));


ang.filter(gtwin(.1,ang.sampleRate));


%% Section 2. Stupid RevCor stuff
myRes = Trial.spk(7);
rang = GetSegs(diff(ang(:,5,7,2)),myRes-61,121,nan);
figure,plot(1:121,nanmean(rang,2))
figure,boundedline(1:121,nanmean(rang,2),2*nanstd(rang,[],2));


%% Section 3 Testing MTAApfs numIter arg

units = [2,7,10,14,20,23,27];
bspfsw = MTAApfs(Trial,units,'walk',overwrite,'spkShuffle',0.5,'numIter',1000);
