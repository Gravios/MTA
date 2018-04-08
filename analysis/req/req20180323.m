

% Spike characterization
% Spike Constellation

dbstop in select_units at 15

select_units([],[],'set');

figure,
subplot(121);  plot([nq.SpkWidthC],[nq.SpkWidthR],'.');
subplot(122);  plot([nq.SpkWidthL],[nq.SpkWidthR],'.');

figure,
subplot(121);  plot([nq.RightMax],[nq.SpkWidthR],'.');
subplot(122);  plot([nq.RightMax],[nq.SpkWidthC],'.');

figure
subplot(121);  plot([nq.SpkWidthR],[nq.CenterMax],'.');
subplot(122);  plot([nq.SpkWidthR],[nq.RightMax],'.');

figure
subplot(121);  plot([nq.SpkWidthR],[nq.AmpSym],'.');
subplot(122);  plot([nq.TimeSym],[nq.SpatLocal],'.');


figure
subplot(151);  plot([nq.TimeSym],[nq.SpkWidthL],'.');
subplot(152);  plot([nq.TimeSym],[nq.SpkWidthC],'.');
subplot(153);  plot([nq.TimeSym],[nq.SpkWidthR],'.');
subplot(154);  plot([nq.TimeSym],log10([nq.FirRate]),'.');
subplot(155);  plot(log10([nq.FirRate]),log10([nq.troughSD]),'.');

figure
plot3([nq.TimeSym],[nq.SpkWidthR],[nq.SpkWidthC],'.');



fet = nunity([[nq.TimeSym],[nq.SpkWidthR],[nq.SpkWidthC],[nq.FirRate],[nq.AmpSym]]);
nind = nniz(fet);

map = tsne(fet(nind,:),[],2,2,80);

figure,plot(map(:,1),map(:,2),'.');

figure,plot(fet(:,1),fet(:,end),'.');

figure();  plot(nq.AmpSym(nind),map(:,1),'.');
figure();  plot(nq.TimeSym(nind),map(:,1),'.');
figure();  plot(nq.FirRate(nind),map(:,1),'.');
figure();  plot(nq.SpkWidthR(nind),map(:,2),'.');
figure();  plot(nq.SpkWidthR(nind),nq.FirRate(nind),'.');
figure();  plot(nq.SpkWidthR(nind),nq.AmpSym(nind),'.');

figure,plot(nq.eDist,nq.SNR,'.');




% Spike constellation

Session = MTATrial.validate('jg05-20120317.cof.all');

lfp = Session.load('lfp',1:96);
tlfp = lfp.copy();
lfp = tlfp.copy();
lfp.filter('ButFilter',3,30,'low');
Session.load('stc','msnn_ppsvd_raux');

stc = label_ripples(Session,[],49:56);

rper = Session.stc{'g'};
rper = rper+[-0.01,0.01];

figure,plot(lfp(rper(4,:),57:64))

groups = {};
for g = 1:numel(Session.parameters.anatomicalDescription.channelGroups.group)
    groups{g} = cell2mat(cf(@(c) str2num(c.Text)+1,...
                         Session.parameters.anatomicalDescription.channelGroups.group{g}.channel));
end

g = 7;
tlfp.data = lfp(:,groups{g});
slfp = permute(tlfp.segs(rper(:,1)-10,60),[2,1,3]);
dslfp = sq(mean(sum(diff(slfp,1,2),2)));
figure,plot(dslfp)




spk = Session.spk.copy();
spk.load_spk(Session);
gid = 7;
spkg = Session.spk.copy();
spkg.clu = spk.clu(ismember(spk.clu,spk.map(spk.map(:,2)==gid,1)));
spkg.res = spk.res(ismember(spk.clu,spk.map(spk.map(:,2)==gid,1)));
spkg.spk = spk.spk(ismember(spk.clu,spk.map(spk.map(:,2)==gid,1)),1:numel(groups{gid}),:);
spkp = sqrt(sum(spk.spk(ismember(spk.clu,spk.map(spk.map(:,2)==gid,1)),1:numel(groups{gid}),:).^2,3));
spkg.sampleRate = spk.sampleRate;

gunits = unique(spkg.clu);

figure,plot(sq(nanmean(spkg.spk(10,:,:))))

figure,hold on,
for g = gunits'
plot(mean(spkg.spk(gunits(g)==spkg.clu,:)))
end


figure,hold on,
for g = gunits',
    plot(ppval(pchip(1:8,mean(spkp(gunits(g)==spkg.clu,1:8))),linspace(1,8,64)))
end

figure,plot(spline(1:8,mean(spkp(gunits(g)==spkg.clu,:))))

figure,hold on,
mspkp = 2e4;%prctile(spkp(:),99);
for g = 1:numel(gunits)
    clup = mean(spkp(gunits(g)==spkg.clu,:));
    vp(g) = (sum(clup(1:4))-sum(clup(5:8)))./sum(clup);
    hp(g) = (sum(clup(1:2:8))-sum(clup(2:2:8)))./sum(clup);
    zp(g) = (1-sum(clup)/mspkp);
    plot3(hp(g),vp(g),zp(g),'.')
end



figure,
for g = 1:numel(gunits)
    clf();
subplot(121),hold on,
s = [];
o = 1;
s(:,o) = ppval(pchip(o:2:8,mean(spkp(gunits(g)==spkg.clu,o:2:8))),linspace(1,8,64));
plot(s(:,o));
o = 2;
s(:,o) = ppval(pchip(o:2:8,mean(spkp(gunits(g)==spkg.clu,o:2:8))),linspace(1,8,64));
plot(s(:,o));
plot(mean(s,2));
subplot(122);
plot(bsxfun(@plus,sq(mean(spks(gunits(g)==spkg.clu,:,:))),fliplr(linspace(1,1e4,8))')')
waitforbuttonpress();
end




% multi session spk constellations
electrodeSet = 'Buz64_1';

n = 1
probes(n).name        = 'Buz64';
probes(n).numChannels = 64;
probes(n).elcGroups   = logical([numElcGroups,probes(n).numChannels]);
probes(n).spkGroups   = logical([numSpkGroups,probes(n).numChannels]);

n = 2
probes(n).name        = 'Lin32';
probes(n).numChannels = 32;
probes(n).elcGroups   = logical([numElcGroups,probes(n).numChannels]);
probes(n).spkGroups   = logical([numSpkGroups,probes(n).numChannels]);

numGroups = size(probes(n).elcGroups(m),1);

Sessions = {MTATrial.validate('jg05-20120316.cof.all'),MTATrial.validate('jg05-20120317.cof.all')};
           cf(@(s)  s.update_parameters(),            Sessions);
           cf(@(s)  s.load('stc','msnn_ppsvd_raux'),  Sessions);
           cf(@(s)  label_ripples(s,[],49:56),        Sessions);
lfp      = cf(@(s)  s.load('lfp',49:56),              Sessions);
           cf(@(l)  l.filter('ButFilter',3,30,'low'), lfp);           
rper     = cf(@(s)  s.stc{'g'}+[-0.01,0.01],          Sessions);
pyrz     = cf(@(l,p) sq(mean(sum(diff(permute(l.segs(p(:,1)-10,60),[2,1,3]),1,2),2))), lfp,rper);

figure();
hold('on');
for i = 1:numel(pyrz)
    plot(pyrz{i});
end

spk      = cf(@(s)    s.spk.copy(),                    Sessions);
spk      = cf(@(k,s)  k.load_spk(s),                   spk,Sessions);

units    = cf(@(k)    k.map(k.map(:,2)==6,1),          spk);


spkg     = cf(@(s)    s.spk.copy(),                    Sessions);
           cf(@(k,s,u)  set(k,'clu',s.clu(ismember(s.clu,u))), spkg,spk,units);
           cf(@(k,s,u)  set(k,'res',s.res(ismember(s.clu,u))), spkg,spk,units);
           cf(@(k,s,u)  set(k,'spk',s.spk(ismember(s.clu,u),1:8,:)), spkg,spk,units);
           cf(@(k,s)  set(k,'sampleRate',s.sampleRate), spkg,spk);
spkp     = cf(@(k)    sqrt(sum(k.spk.^2,3)), spkg);



figure,hold on,
uclr = 'rb';
mspkp = 2e4;%prctile(spkp(:),99);
for s = 1:2,
    for g = 1:numel(units{s})
        clup = mean(spkp{s}(units{s}(g)==spkg{s}.clu,:));
        vp(g) = (sum(clup(1:4))-sum(clup(5:8)))./sum(clup);
        hp(g) = (sum(clup(1:2:8))-sum(clup(2:2:8)))./sum(clup);
        zp(g) = (1-sum(clup)/mspkp);
        %plot3(hp(g),vp(g),zp(g),'.b')
        plot(zp(g),vp(g),['.',uclr(s)])    
    end
end




figure,
for g = 1:numel(gunits)
    clf();
subplot(121),hold on,
s = [];
o = 1;
s(:,o) = ppval(pchip(o:2:8,mean(spkp(gunits(g)==spkg.clu,o:2:8))),linspace(1,8,64));
plot(s(:,o));
o = 2;
s(:,o) = ppval(pchip(o:2:8,mean(spkp(gunits(g)==spkg.clu,o:2:8))),linspace(1,8,64));
plot(s(:,o));
plot(mean(s,2));
subplot(122);
plot(bsxfun(@plus,sq(mean(spks(gunits(g)==spkg.clu,:,:))),fliplr(linspace(1,1e4,8))')')
waitforbuttonpress();
end



Sessions = {MTATrial.validate('jg05-20120316.cof.all'),MTATrial.validate('jg05-20120317.cof.all')};
for t = 1:numel(Sessions),
    s = MTASession.validate(Sessions{t});
    s.spk.create(s);
    s.save();
end
Sessions = {MTATrial.validate('jg05-20120316.cof.all'),MTATrial.validate('jg05-20120317.cof.all')};
spk      = cf(@(s)    s.spk.copy(),                    Sessions);
spk      = cf(@(k,s)  k.load_spk(s),                   spk,Sessions);

units    = cf(@(k)    k.map(k.map(:,2)==6,1),          spk);

pfs = cf(@(s,u) pfs_2d_states(s,u,[],[],[],[],true,1),  Sessions, units);


nx = numel(pfs{1});
ny = numel(pfs);



spk{1}.map(spk{1}.map(:,2)==8,:)
spk{2}.map(spk{2}.map(:,2)==8,:)




% Testing Pair


unit = {14,33};

hfig = figure(39394);
set(hfig,'Position',[       97         528        1460         207])
clf();
for y = 1:ny,
    for x = 1:nx,
        subplot2(ny,nx,y,x);
        plot(pfs{y}{x},unit{y},'mean',true,[],true);
    end
end


    
    
    
    
    


