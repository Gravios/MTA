MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
Trial = MTATrial('jg05-20120317');
fwin = gausswin(11)./sum(gausswin(11));
Trial.load('ang');
Trial.load('xyz');
Trial.xyz.filter(fwin);

vel =[];
vel = [vel,Filter0(gausswin(61)./sum(gausswin(61)),log10(Trial.vel({'spine_lower'}))).^2];
vel = [vel,Filter0(fwin,log10(Trial.vel({'spine_lower','spine_upper','head_front'},[1,2])))];
vel = [vel,Filter0(fwin,log10(Trial.vel({'spine_lower','spine_upper','head_front'},3)))];
vel = [vel,Filter0(gausswin(61)./sum(gausswin(61)),sqrt(sum(diff(Trial.com(Trial.xyz.model.rb({'spine_lower','pelvis_root','spine_middle'}))).^2,3)))];  
vel = clip(vel,-2,2);
vel = [vel,Trial.ang(1:end-1,1,2,2)];
vel = [vel,Trial.ang(1:end-1,3,4,2)];
vel = [vel,vel.^2];
vel(isnan(vel)) = 0;

[State, hmm, decode] = gausshmm(vel,12); % Failed kinda

Trial.ang.filter(gausswin(11)./sum(gausswin(11)));

ind = Trial.xyz(:,1,3)~=0;
fpv = cat(1,[-2,-2],log10(Filter0(gausswin(61)./sum(gausswin(61)),Trial.vel({'spine_lower','head_front'}))));
vel =[];
vel = [vel,fpv(ind)];
vel = [vel,Trial.ang(ind,1,2,2)];
vel = [vel,Trial.ang(ind,3,4,2)];
vel = [vel,Trial.ang(ind,1,6,3)];
vel = [vel,Trial.ang(ind,1,8,3)];
vel = [vel,Trial.xyz(ind,7,3)];
vel(isnan(vel)) = 0;

%vel = MTADxyz([],[],Filter0(fwin,vel),Trial.xyz.sampleRate);



[State2, hmm2, decode2] = gausshmm(vel,8);

Stateall = zeros(Trial.xyz.size(1),1);
Stateall(ind) = State2;
figure,plot(Stateall)
Lines(Trial.stc{'r'}.data(:),[],'r');

st = MTADxyz([],[],Stateall',Trial.xyz.sampleRate);



name = 'hmm_test20140311';

Stc = MTAStateCollection(Trial.spath,Trial.filebase,name,Trial.stc.sync.copy,Trial.stc.origin,1,'stc');




g = 1;l = 2;
Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   fasts==g,...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'hwalk','g','TimeSeries');
Stc{'g'}.cast('TimePeriods');

rs = zeros(st.size(1),1);
rs(ismember(st.data,[3,7,8])) = 1;