function bhv_rhm_spk(Trial,varargin)
[stc_mode,autoIncr,overwrite] = DefaultArgs(varargin,{'auto_wbhr',false,false});

%sname = 'jg04-20120129';
%sname = 'jg05-20120317';
%sname = 'jg05-20120309';
%sname = 'er01-20110719';
%tname = 'all';
%stc_mode = 'auto_wbhr';
%Trial = MTATrial(sname,tname);


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

try load(fullfile(Trial.path.cfg,'super_fet_rhm.normparm.mat'));end
if ~exist('rhm_mean','var'),
rhm_std  = nanstd(rhm(vnn,:));
rhm_mean = nanmean(rhm(vnn,:));
save(fullfile(Trial.path.cfg,'super_fet_rhm.normparm.mat'),'rhm_mean','rhm_std');
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


%sind = Trial.stc{'w',ar.sampleRate};

%overwrite = true;

tind = Trial.stc{'t',ar.sampleRate};
units = select_units(Trial,25,'pyr');


pfar = MTAApfs(Trial,units,tind,overwrite,[],[.05,.05],[1.8,1.8],'type','ar','xyz',ar,'bound_lims',bound_limx);
pfav = MTAApfs(Trial,units,tind,overwrite,[],[.05,.05],[1.8,1.8],'type','av','xyz',av,'bound_lims',bound_limv);
pfrv = MTAApfs(Trial,units,tind,overwrite,[],[.05,.05],[1.8,1.8],'type','rv','xyz',rv,'bound_lims',bound_limr);

%dthresh = mean(pfx.parameters.binDims./2*sqrt(2));
dthresh = 30;
pfarsx={};
pfavsx={};
pfrvsx={};
stss='twr';
uSampInd = [];

for s = 1:numel(stss);
    sind = Trial.stc{stss(s),ar.sampleRate};
    pfx{s} = MTAApfs(Trial,units,sind,overwrite);
    for i = units
        halfMaxRate = max(pfx{1}.data.rateMap(:,i==pfx{1}.data.clu))/4;
        sxind = find(pfx{s}.data.rateMap(:,i==pfx{s}.data.clu)>halfMaxRate);
        sxsubs = Ind2Sub(pfx{s}.adata.binSizes',sxind);
        sxpos = [pfx{s}.adata.bins{1}(sxsubs(:,1)),pfx{s}.adata.bins{2}(sxsubs(:,2))];
        sxind = false([xyz.size(1),1]);
        for j = 0:(prod(pfx{1}.adata.binSizes)/10);
            if (j*10+10)>size(sxpos,1),
                eind = size(sxpos,1)-j*10;
            else
                eind = 10;
            end
            sxind = sxind|sum(sqrt(sum((repmat(permute(sxpos((j*10+1):(j*10+eind),:),[3,1,2]),[xyz.size(1),1,1])-repmat(xyz(:,Trial.trackingMarker,[1,2]),[1,eind,1])).^2,3))<dthresh,2)>=1;
        end
        if s==1,uSampInd(:,i) = sxind;end
        %plot(xyz(sxind,7,1),xyz(sxind,7,2),'.'),xlim([-500,500]),ylim([-500,500])
        %hold on
        %plot(sxpos(:,1),sxpos(:,2),'.r'),xlim([-500,500]),ylim([-500,500])
        if sum(sxind)==0,
            pfarsx{i,s}={};
            pfavsx{i,s}={};
            pfrvsx{i,s}={};
            continue,
        else
            nind = sind&ThreshCross(sxind,.5,3);        
        end

        pfarsx{i,s} = MTAApfs(Trial,i,nind,overwrite,[stss(s) 'mrsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                              'type','ar','xyz',ar,'bound_lims',bound_limx);
        pfavsx{i,s} = MTAApfs(Trial,i,nind,overwrite,[stss(s) 'mvsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                              'type','av','xyz',av,'bound_lims',bound_limv);
        pfrvsx{i,s} = MTAApfs(Trial,i,nind,overwrite,[stss(s) 'mrvsUnit' num2str(i)],[.05,.05],[1.8,1.8],...
                              'type','rv','xyz',rv,'bound_lims',bound_limr);
        
    end
end



[accg,tbin] = autoccg(Trial,[],'theta');


ny =4;
hFig = figure(858283);
set(hFig,'position',[165,31,913,740]);
i=units(1);
while i~=-1,
    subplot2(ny,numel(stss)+1,1,1);
    try,plot(xyz(uSampInd(:,i),7,1),xyz(uSampInd(:,i),7,2),'.'),xlim([-500,500]),ylim([-500,500]),end
    subplot2(ny,numel(stss)+1,2,1);
    bar(tbin,accg(:,i)),axis tight,
    for s = 1:numel(stss)
        try,
            subplot2(ny,numel(stss)+1,1,s+1);
            pfx{s}.plot(i,[],1);
            title(num2str(i));
            set(gca,'YTickMode','manual');
            set(gca,'yticklabel',{});
            set(gca,'xTickMode','manual');
            set(gca,'xticklabel',{});
        end
        try,
            subplot2(ny,numel(stss)+1,2,s+1);
            pfrvsx{i,s}.plot(i,[],1);
            text(bound_limr(1,1)+.2,bound_limr(2,1)+.2,sprintf('%2.1f',nanmax(pfrvsx{i,s}.data.rateMap(:,i==pfrvsx{i,s}.data.clu))),'Color','w','FontWeight','bold','FontSize',10);
            hTxt = text(bound_limr(1,1)+.2,bound_limr(2,1)+1,'log10(vel)','Color','w','FontWeight','bold','FontSize',8);
            set(hTxt, 'rotation', 90)
            text(bound_limr(1,1)+1,bound_limr(2,1)+.2,'rhm(6-12)pow','Color','w','FontWeight','bold','FontSize',8);
        end
        try,
            subplot2(ny,numel(stss)+1,3,s+1);
            pfarsx{i,s}.plot(i,[],1);            
           text(bound_limx(1,1)+.2,bound_limx(2,1)+.2,sprintf('%2.1f',nanmax(pfarsx{i,s}.data.rateMap(:,i==pfarsx{i,s}.data.clu))),'Color','w','FontWeight','bold','FontSize',10);
            hTxt = text(bound_limx(1,1)+.2,bound_limx(2,1)+1,'rhm(6-12)pow','Color','w','FontWeight','bold','FontSize',8);
            set(hTxt, 'rotation', 90)
            text(bound_limx(1,1)+1,bound_limx(2,1)+.2,'angle radians','Color','w','FontWeight','bold','FontSize',8);
        end
        try,
            subplot2(ny,numel(stss)+1,4,s+1);
            pfavsx{i,s}.plot(i,[],1);
            text(bound_limv(1,1)+.2,bound_limv(2,1)+.2,sprintf('%2.1f',nanmax(pfavsx{i,s}.data.rateMap(:,i==pfavsx{i,s}.data.clu))),'Color','w','FontWeight','bold','FontSize',10);
            hTxt = text(bound_limv(1,1)+.2,bound_limv(2,1)+1,'log10(vel)','Color','w','FontWeight','bold','FontSize',8);
            set(hTxt, 'rotation', 90)
            text(bound_limv(1,1)+1,bound_limv(2,1)+.2,'angle radians','Color','w','FontWeight','bold','FontSize',8);
        end
    end
    i = figure_controls(hFig,i,units,autoIncr);
    pause(.1)
    reportfig(fullfile(Trial.path.data,'figures'),hFig, ...
          ['BHV_SPK_ang_vel_xy'],[],[Trial.filebase ' :unit:' num2str(i)],200,true);

end






