addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/fwdrebayesiandecodingtools/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/tempScripts/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/aux/

MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;
trialList = 'Ed10VR_teleport';
OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';

T = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');

Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);

Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

overwrite = false;
units =[1,5,7,9,16,18,21,22,23,24,25,28,29,31,38,39,43,49,51,52,73,74,77,79,80,85,88,94,99,101,104,134,138,142,158,168,184,185];

% Generate unit auto correlogram
[accg,tbin] = autoccg(Trial,units,'theta');


binDims = [40,40];
smoothingWeights = [1.2,1.2];
pfs = {};
for t = 1:nt-1
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial); 
   for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite, ...
                           'binDims',binDims,'SmoothingWeights',smoothingWeights);
    end
end
t = nt;
Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
Trial.stc = Stc.copy;
Trial.stc.load(Trial); 
for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,[states{i},'-shifted'],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'type','xy');
end
t = nt+1;
for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,['shifted&',states{i}],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'type','xy');
end


units =[1,5,7,9,16,18,21,22,23,24,25,28,29,31,38,39,43,49,51,52,73,74,77,79,80,85,88,94,99,101,104,134,138,142,158,168,184,185];
spCohere = [];
for u = 1:numel(units),
    spCohere(end+1) = pfs{2,1}.spatialCoherence(units(u));
end


sunits = units(spCohere>0.985);


ratemap=[];
for u = 1:numel(sunits),
    ratemap(:,:,u) = pfs{2,1}.plot(sunits(u),'isCircular',false);
end






Trial = MTATrial.validate(T(1));
Trial.stc = Stc.copy;
Trial.stc.load(Trial); 
xyz = Trial.load('xyz');
xyz.resample(10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'theta',sunits,1.5,true);

E = decode_bayesian_poisson(ratemap,ufr.data');

xbins = pfs{2,1}.adata.bins{1};
ybins = pfs{2,1}.adata.bins{2};

pos = nan([size(E,3),2]);
for tind = 1:size(E,3),
    try,
        xyi = LocalMinimaN(-E(:,:,tind),0,100);        
        pos(tind,:) = [xbins(xyi(1)),ybins(xyi(2))];
    end
end    

figure,plot(sq(xyz(:,5,[1])))
hold on,
plot(pos(:,1).*double(sum(ufr.data,2)>100));


xp = sq(xyz(:,5,[1,2]))-pos;
[th,r] = cart2pol(xp(:,1),xp(:,2));



figure,plot(r(sum(ufr.data,2)>1))


figure,
plot(r)
hold on
plot(sum(ufr.data,2))


figure,plot(sq(xyz(:,5,[1])))
hold on,
[~,ttt] = min(-mean(E,2));
%ttt(ttt==0)=1;
plot(xbins(sq(ttt)));



xd = sq(xyz(:,5,[1]))-pos(:,1);

figure,
plot(xd.*double(sum(ufr.data,2)>10))
hold on
plot(sum(ufr.data,2))


tind = 800:1300;
figure,

plot(xyz(tind,5,1),xyz(tind,5,2),'.')
hold on
plot(pos(tind,1).*sum(ufr.data(tind,:),2)~=0),pos(tind,2).*sum(ufr.data(tind,:),2)==0),'.')
