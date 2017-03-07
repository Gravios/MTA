function bhv_ncpp_JPDF(Trial,varargin)
[mode,stc_mode,marker] = DefaultArgs(varargin,{{'height','hangle'},'auto_wbhr','spine_lower'});

% $$$ sname = 'Ed05-20140528'
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


figH = figure(22030233); 
set(figH,'pos',[14,325,1181,420+(420*(numel(mode)-1))]);



[ncp,fs] = fet_ncp(Trial,[],'Swspectral');
ncp.data  = log10(ncp.data);
ncp.data(~nniz(ncp.data(:)))=nan;

vel = xyz.vel(1);
vel.resample(ncp);
vnn = nniz(vel);

% $$$ try load(fullfile(Trial.path.cfg,'super_fet_ncp.normparm.mat'));end
% $$$ if ~exist('ncp_mean','var'),
ncp_std  = nanstd(ncp(vnn,:));
ncp_mean = nanmean(ncp(vnn,:));
% $$$ save(fullfile(Trial.path.cfg,'super_fet_ncp.normparm.mat'),'ncp_mean','ncp_std');
% $$$ end

ncp.data = (ncp.data-repmat(ncp_mean,[ncp.size(1),1]))./repmat(ncp_std,[ncp.size(1),1]);


ncp_snr = nanmedian(ncp(:,fs>5&fs<13),2).*abs(nanmedian(ncp(:,fs>5&fs<13),2)-nanmedian([ncp(:,fs<4),ncp(:,fs>14&fs<18)],2));
ncp_snr(ncp_snr<0)=.001;
ncp_snr = log10(ncp_snr);
ncp.data = ncp_snr;


for m = 1:numel(mode),

switch mode{m}
  case 'height'
    vel = MTADxyz('data',xyz(:,'head_front',3),'sampleRate',xyz.sampleRate);
    vel.resample(ncp);
    vel.data = log10(vel.data);
    vnn = nniz(vel);
    xlab = 'Log10 Head Height (mm)';
    vlim =[1.2,2.4];
  case 'hangle'
    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    vel = MTADxyz('data',ang(:,'head_back','head_front',2),'sampleRate',xyz.sampleRate);
    vel.resample(ncp);
    vnn = nniz(vel);
    xlab = 'head angle (rad)';
    vlim =[-1.5,1.5];
  case 'bspeed'
    vel = xyz.vel(marker,[1,2]);
    vel.resample(ncp);
    vel.data = log10(vel.data);
    vnn = nniz(vel);
    xlab = 'Log10 Body Speed (cm/s)';
    vlim =[-.5,1.5];
end



s = 'w';
sncp = ncp(Trial.stc{s},:);
svel = vel(Trial.stc{s},:);

vbins = 50;
vedgs = linspace(vlim(1),vlim(2),vbins);
[N,vbs] = histc(svel,vedgs);

subplot2(numel(mode),2,m,1);
hist2([svel,sncp],linspace(vlim(1),vlim(2),50),linspace(-2,1.3,70)),
title(['NCP ' mode{m} ' JPDF']);
xlabel(xlab);
ylabel('NCP super fet');

subplot2(numel(mode),2,m,2);
bar(vedgs,N,'histc')
title(['Marginal Distrb of ' mode{m}]);
xlabel(xlab);
ylabel('count');

end

suptitle([Trial.filebase ' : ' Trial.stc{s}.label]);
reportfig(fullfile(Trial.path.data,'figures'),figH, ...
          ['NCP_' strjoin(mode,'_') '_JPDF' ],[],[Trial.filebase ' :' Trial.stc{s}.label],200,true);

