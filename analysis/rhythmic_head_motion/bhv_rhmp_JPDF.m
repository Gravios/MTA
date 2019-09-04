function bhv_rhmp_JPDF(Trial,varargin)
[mode,stc_mode,marker] = DefaultArgs(varargin,{{'height','hangle'},'auto_wbhr','spine_lower'});

% $$$ sname = 'jg05-20120311';
% $$$ tname = 'all';
% $$$ mode = {'height','hangle'};
% $$$ marker = 'spine_lower';
% $$$ stc_mode = 'auto_wbhr';
% $$$ Trial = MTATrial(sname,tname);

Trial.stc.updateMode(stc_mode);
Trial.stc.load;


xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.2,xyz.sampleRate));


figH = figure(22030232); 
set(figH,'pos',[14,325,1181,420+(420*(numel(mode)-1))]);



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


rhm_snr = nanmedian(rhm(:,fs>6&fs<12),2).*abs(nanmedian(rhm(:,fs>6&fs<12),2)-nanmedian([rhm(:,fs<4),rhm(:,fs>14&fs<16)],2));
rhm_snr(rhm_snr<0)=.001;
rhm_snr = log10(rhm_snr);
rhm.data = rhm_snr;


for m = 1:numel(mode),

switch mode{m}
  case 'height'
    vel = MTADxyz('data',xyz(:,'head_front',3),'sampleRate',xyz.sampleRate);
    vel.resample(rhm);
    vel.data = log10(vel.data);
    vnn = nniz(vel);
    xlab = 'Log10 Head Height (mm)';
    vlim =[1.2,2.4];
  case 'hangle'
    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    vel = MTADxyz('data',ang(:,'head_back','head_front',2),'sampleRate',xyz.sampleRate);
    vel.resample(rhm);
    vnn = nniz(vel);
    xlab = 'head angle (rad)';
    vlim =[-1.5,1.5];
  case 'bspeed'
    vel = xyz.vel(marker,[1,2]);
    vel.resample(rhm);
    vel.data = log10(vel.data);
    vnn = nniz(vel);
    xlab = 'Log10 Body Speed (cm/s)';
    vlim =[-.5,1.5];
end



s = 'w';
srhm = rhm(Trial.stc{s},:);
svel = vel(Trial.stc{s},:);

vbins = 50;
vedgs = linspace(vlim(1),vlim(2),vbins);
[N,vbs] = histc(svel,vedgs);

subplot2(numel(mode),2,m,1);
hist2([svel,srhm],linspace(vlim(1),vlim(2),50),linspace(-2,1,50)),
title(['RHM ' mode{m} ' JPDF']);
xlabel(xlab);
ylabel('RHM super fet');

subplot2(numel(mode),2,m,2);
bar(vedgs,N,'histc')
title(['Marginal Distrb of ' mode{m}]);
xlabel(xlab);
ylabel('count');

end

suptitle([Trial.filebase ' : ' Trial.stc{s}.label]);
reportfig(fullfile(Trial.path.project,'figures'),figH, ...
          ['RHM_' strjoin(mode,'_') '_JPDF' ],[],[Trial.filebase ' :' Trial.stc{s}.label],200,true);

