sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

t = 2;
Trial = Trials{t};
xyz = Trial.load('xyz');
stc = Trial.stc.copy;

spk = Trial.spk.copy;
spk.create(Trial,xyz.sampleRate,'loc',[],'deburst');

spkt = Trial.spk.copy;
spkt.create(Trial,xyz.sampleRate,'theta-groom-sit',[],'deburst');

% GET theta phase
lfp = Trial.load('lfp',sessionList(t).thetaRef);
lfp.resample(xyz);
thetaPhase = lfp.phase([5,13]);

% GET placefield stats and analysis subset
% $$$ pfstats = compute_pfstats_bs(Trial,'overwrite',false);
% $$$ units = pfstats.cluMap;

ang = create(MTADang,Trial,xyz);

pft = pfs_2d_theta(Trial,'overwrite',true);

%pft.data.si(86)

units = select_placefields(Trial);
sWeights = 1:0.5:4;
filtCutOffFreq = [0.5,1,1.5,2.5,5,10,20];


pf = {};  drz = {};  ddz = {};
for s = 1:numel(sWeights),
    defargs = get_default_args('MjgER2016','MTAApfs','struct');
    defargs.units = units;
    defargs.states = 'theta-groom-sit';
    defargs.SmoothingWeights = repmat(sWeights(s),[1,2]);
    defargs.overwrite = overwrite;
    defargs = struct2varargin(defargs);
    pf{s} = MTAApfs(Trial,defargs{:});
    for f = 1:numel(filtCutOffFreq),
        drz{s,f} = compute_drz(Trial,pf{s},units,[],filtCutOffFreq(f));    
        %ddz{s,f} = compute_ddz(Trial,pf{s},units,[],filtCutOffFreq(f));
    end
end


unit = 86;

res = spk(unit);
%res = spkt(unit);
res(res>xyz.size(1))=[];
phzspk = thetaPhase(res,spk.map(spk.map(:,1)==unit,2));

figure(2017120501);clf();
nx = numel(sWeights);
ny = numel(filtCutOffFreq)+1;
for s = 1:numel(sWeights),
    subplot2(ny,nx,1,s);  plot(pf{s},unit);    
    for f = 1:numel(filtCutOffFreq),
        drzspk = drz{s,f}(res,unit==units);
        gind = ~isnan(drzspk)&~isnan(phzspk);                
        
        subplot2(ny,nx,f+1,s);
        if sum(gind)>10,
            hold('on');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
            plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
            xlim([-1,1]),
            ylim([-180,540])
        end                
    end
end


s = 5;f = 7; 
figure();
drzspk = drz{s,f}(res,unit==units);
gind = ~isnan(drzspk)&~isnan(phzspk);                
hold('on');
plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
xlim([-1,1]),
ylim([-180,540])

nbins = 100;
csi = hsv(nbins);
csi = discretize(ang(res,5,7,1),linspace(-pi,pi,nbins));
cs = jet(nbins);
csi = discretize(ang(res,5,7,2),linspace(-pi/2,pi/2,nbins));

figure;hold on
scatter(drzspk(gind),circ_rad2ang(phzspk(gind)),50,cs(csi(gind),:),'filled');
scatter(drzspk(gind),circ_rad2ang(phzspk(gind))+360,50,cs(csi(gind),:),'filled');


figure(2017120502);clf();
nx = numel(sWeights);
ny = numel(filtCutOffFreq)+1;
for s = 1:numel(sWeights),
    subplot2(ny,nx,1,s);  plot(pf{s},unit);    
    for f = 1:numel(filtCutOffFreq),
        ddzspk = ddz{s,f}(res,unit==units);        
        gind = ~isnan(ddzspk)&~isnan(phzspk);
        subplot2(ny,nx,f+1,s);
        if sum(gind)>10,
            hold('on');            
            plot(ddzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
            plot(ddzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
            xlim([-300,300]),
            ylim([-180,540])
        end
    end
end


figure();
s = 4;
f = 4;
ddzspk = ddz{s,f}(res,unit==units);
gind = ~isnan(ddzspk)&~isnan(phzspk);
hold('on');            
plot(ddzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
plot(ddzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
xlim([-300,300]),
ylim([-180,540])
