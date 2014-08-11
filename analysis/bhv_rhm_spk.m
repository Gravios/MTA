
%Psname = 'jg04-20120129';
sname = 'jg05-20120317';
sname = 'jg05-20120310';
%sname = 'er01-20110719';
tname = 'all';
mode = {'height','hangle'};
marker = 'spine_lower';
stc_mode = 'auto_wbhr';
Trial = MTATrial(sname,tname);


Trial.stc.updateMode(stc_mode);
Trial.stc.load;
Trial.stc.states{5} = theta(Trial);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.2,xyz.sampleRate));


[rhm,fs] = fet_rhm(Trial,xyz.sampleRate,'Swspectral');
rhm.data  = log10(rhm.data);
rhm.data(rhm<-8) = nan;
rhm.data(~nniz(rhm.data(:)))=nan;
vel = xyz.vel(1);
vel.resample(rhm);
vnn = nniz(vel);

try load(fullfile(Trial.path.MTAPath,'super_fet_rhm.normparm.mat'));end
if ~exist('rhm_mean','var'),
rhm_std  = nanstd(rhm(vnn,:));
rhm_mean = nanmean(rhm(vnn,:));
save(fullfile(Trial.path.MTAPath,'super_fet_rhm.normparm.mat'),'rhm_mean','rhm_std');
end
rhm.data = (rhm.data-repmat(rhm_mean,[rhm.size(1),1]))./repmat(rhm_std,[rhm.size(1),1]);

% $$$ rhm_snr = nanmedian(rhm(:,fs>6&fs<12),2).*abs(nanmedian(rhm(:,fs>6&fs<12),2)-nanmedian([rhm(:,fs<4),rhm(:,fs>14&fs<16)],2));
% $$$ rhm_snr(rhm_snr<0)=.001;
% $$$ rhm_snr = log10(rhm_snr);
% $$$ rhm.data = rhm_snr;
rhm.data = nanmedian(rhm(:,fs>6&fs<12),2);

ufr = Trial.ufr.copy;
ufr.create(Trial,xyz,'twin',.01);
ufr.filter(gtwin(.4,xyz.sampleRate));
ufr.resample(rhm);

ang = Trial.ang.copy;
ang.create(Trial,xyz);
ang.data = ang.data(:,5,7,2);
ang.resample(rhm);

xyz.resample(rhm);

vel = xyz.vel(7);
bound_limx = [-1.4,1.4;-1,1.75];
ar = MTADxyz('data',[ang.data,rhm.data],'sampleRate',rhm.sampleRate);

bound_limv = [-1.4,1.4;-2,2];
av = MTADxyz('data',[ang.data,log10(vel.data)],'sampleRate',rhm.sampleRate);

bound_limr = [-1,1.75;-2,2];
rv = MTADxyz('data',[rhm.data,log10(vel.data)],'sampleRate',rhm.sampleRate);

sind = Trial.stc{'w',ar.sampleRate};
%sind = Trial.stc{'t',ar.sampleRate};


units = select_units(Trial,25,'pyr');

pfar = MTAApfs(Trial,units,sind,true,[],[.05,.05],[1.8,1.8],'type','ar','xyz',ar,'bound_lims',bound_limx);
pfav = MTAApfs(Trial,units,sind,true,[],[.05,.05],[1.8,1.8],'type','av','xyz',av,'bound_lims',bound_limv);
pfrv = MTAApfs(Trial,units,sind,true,[],[.05,.05],[1.8,1.8],'type','rv','xyz',rv,'bound_lims',bound_limr);

pfx = MTAApfs(Trial,units,sind,true);

%dthresh = mean(pfx.parameters.binDims./2*sqrt(2));
dthresh = 30;
pfarsx={};
pfavsx={};
for i = units
    halfMaxRate = max(pfx.data.rateMap(:,i==pfx.data.clu))/4;
    sxind = find(pfx.data.rateMap(:,i==pfx.data.clu)>halfMaxRate);
    sxsubs = Ind2Sub(pfx.adata.binSizes',sxind);
    sxpos = [pfx.adata.bins{1}(sxsubs(:,1)),pfx.adata.bins{2}(sxsubs(:,2))];
    sxind = false([xyz.size(1),1]);
    for j = 0:(prod(pfx.adata.binSizes)/10);
        if (j*10+10)>size(sxpos,1),
            eind = size(sxpos,1)-j*10;
        else
            eind = 10;
        end
        sxind = sxind|sum(sqrt(sum((repmat(permute(sxpos((j*10+1):(j*10+eind),:),[3,1,2]),[xyz.size(1),1,1])-repmat(xyz(:,Trial.trackingMarker,[1,2]),[1,eind,1])).^2,3))<dthresh,2)>=1;
    end
    %clf
%    plot(xyz(sxind,7,1),xyz(sxind,7,2),'.'),xlim([-500,500]),ylim([-500,500])
%hold on
%plot(sxpos(:,1),sxpos(:,2),'.r'),xlim([-500,500]),ylim([-500,500])


    nind = sind&ThreshCross(sxind,.5,3);
    pfarsx{i} = MTAApfs(Trial,i,nind,true,['mrsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                      'type','ar','xyz',ar,'bound_lims',bound_limx);
    pfavsx{i} = MTAApfs(Trial,i,nind,true,['mvsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                      'type','av','xyz',av,'bound_lims',bound_limv);
    pfrvsx{i} = MTAApfs(Trial,i,nind,true,['mrvsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                      'type','rv','xyz',rv,'bound_lims',bound_limr);
    
end

[accg,tbin] = autoccg(Trial,[],'theta');

for i = units,
    
    subplot2(4,4,[1:3],1);plot(ang(sind),ufr(sind,i),'.'),xlim(bound_limx(1,:))
    subplot2(4,4,[4],1);hist(ang(sind),linspace(bound_limx(1,1),bound_limx(1,2),60)),xlim(bound_limx(1,:))
    %subplot2(4,4,[1:3],2);plot(rhm(sind),ufr(sind,i),'.'),xlim(bound_limx(2,:))
    %subplot2(4,4,[4],2);hist(rhm(sind),linspace(bound_limx(2,1),bound_limx(2,2),60)),xlim(bound_limx(2,:))
    subplot2(4,4,[1:2],2);plot(rhm(sind),ufr(sind,i),'.'),xlim(bound_limx(2,:))
    subplot2(4,4,[3:4],2);pfrvsx{i}.plot(i,[],1);
    subplot2(4,4,[1:2],3);bar(tbin,accg(:,i)),axis tight,%pfar.plot(i,[],1);
    subplot2(4,4,[3:4],3);pfx.plot(i,[],1);
    subplot2(4,4,[1:2],4);pfarsx{i}.plot(i,[],1);
    subplot2(4,4,[3:4],4);pfavsx{i}.plot(i,[],1);
    title(num2str(i))
waitforbuttonpress
end






