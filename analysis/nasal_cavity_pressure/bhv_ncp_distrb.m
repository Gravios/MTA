function bhv_ncp_distrb(Trial,varargin)
% bhv_ncp_distrb(Trial,varargin)
% behavioral nasal cavity pressure distribution
%
% The mean power spectral density of the nasal cavity pressure sensor
%  binned by various behavioral variables.
%
%  varargin:
%
%      state -    string: label or key of behavioral state periods
%                          you whish to include in the analysis
%
%      mode -     string: the name of the variable comparisons 
%                          see the switch statement in the code
%                          to see which options are available.
%
%      stc_mode - string: the MTAStateCollection which contains the
%                          states you require for the analysis
%
%      marker -   string: the marker name, used for general tracking
%
%      ncpChan -  number: the channel number on which the nasal
%                          cavity pressure signal was recorded
%

[state,mode,stc_mode,marker,ncpChan] = DefaultArgs(varargin,{'walk',{'height','hangle'},'auto_wbhr','spine_lower',2});

% $$$ sname = 'jg05-20120317';
% $$$ tname = 'all';
% $$$ mode = {'height','hangle'};
% $$$ marker = 'spine_lower';
% $$$ stc_mode = 'auto_wbhr';
% $$$ Trial = MTATrial(sname,tname);

Trial.stc.updateMode(stc_mode);
Trial.stc.load;



xyz = Trial.load('xyz');
xyz.filter(gtwin(.2,xyz.sampleRate));


figH = figure(22030233); 
set(figH,'pos',[14,325,1181,420+(420*(numel(mode)-1))]);


[rhm,fs] = fet_ncp(Trial,xyz.sampleRate,'wcsd',ncpChan);
rhm.data  = log10(rhm.data);
rhm.data(rhm<-9) = nan;
rhm.data(nniz(rhm.data))=nan;
vel = xyz.vel(1);
vel.resample(rhm);
vnn = nniz(vel);
rhm.data = (rhm.data-repmat(nanmean(rhm(vnn,:)),[rhm.size(1),1]))./repmat(nanstd(rhm(vnn,:)),[rhm.size(1),1]);


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
  case 'DistMazeCenter'
    xyz.addMarker('maze_center',[0,0,1],{},zeros([xyz.size(1),1,xyz.size(3)]));
    xyz.data(:,xyz.model.gmi('head_front'),3) = 0;
    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    vel = MTADxyz('data',ang(:,'maze_center','head_front',3),'sampleRate',xyz.sampleRate);
    vel.resample(rhm);
    vel.data = log10(vel.data);    
    vnn = nniz(vel);
    [out,x] = MakeUniformDistr(vel(vnn),0,2.6);
    vel.data(vnn) = out;
    xlab = {'Log10 Head and Maze Center Distance (mm)','Note: Distribution made Uniform'};
    vlim =[0,2.6];

end



s = state;
srhm = rhm(Trial.stc{s},:);
svel = vel(Trial.stc{s},:);


vbins = 50;
vedgs = linspace(vlim(1),vlim(2),vbins);
[N,vbs] = histc(svel,vedgs);

mrv = nan(numel(N),rhm.size(2));
for f =1:rhm.size(2),
    mrv(:,f) = accumarray(vbs(nniz(vbs)),srhm(nniz(vbs),f),[vbins,1],@nanmean);
end
mrv(N<20,:) = nan;


subplot2(numel(mode),2,m,1);
%figure
imhand = imagescnan({vedgs,fs,mrv'},prctile(mrv(nniz(mrv(:))),[5,95]),false,1,[0,0,0]);axis xy,
set(imhand,'tag', [mfilename,'-',mode{m}])
title(['Mean NCP psd Given ' mode{m}]);
xlabel(xlab);
ylabel('NCP frequency (Hz)');

subplot2(numel(mode),2,m,2);
bar(vedgs,N,'histc')
title(['Marginal Distrb of ' mode{m}]);
xlabel(xlab);
ylabel('NCP frequency (Hz)');

end

suptitle([Trial.filebase ' : ' Trial.stc{s}.label]);
reportfig(Trial,figH,['NCP_psd_distrib_' strjoin(mode,'_')],[],[Trial.filebase ' :' Trial.stc{state}.label],200,true);





% $$$ spowa = MTADlfp('data',nanmean(rhm(:,fs>6&fs<12),2),'sampleRate',Trial.xyz.sampleRate);
% $$$ 
% $$$ figure,hist2([vel(vnn),spowa(vnn)],75,75)
% $$$ title('all');
% $$$ xlabel('height cm/s')
% $$$ ylabel('log10 rhm power 6-12 hz')
% $$$ 
% $$$ 
% $$$ 
% $$$ splim = prctile(spowa(Trial.stc{'w'},:),[2,98]);
% $$$ vlim  = prctile(  vel(Trial.stc{'w'})  ,[2,98]);
% $$$ 
% $$$ s = 'w';
% $$$ vbins = 25;
% $$$ vedgs = linspace(vlim(1),vlim(2),vbins);
% $$$ spbins = 25;
% $$$ spedgs = linspace(splim(1),splim(2),spbins);
% $$$ figure,hist2([vel(Trial.stc{s}),spowa(Trial.stc{s})],vedgs,spedgs)
% $$$ title(Trial.stc{s}.label);
% $$$ xlabel('height cm/s')
% $$$ ylabel('log10 rhm power 6-12 hz')


%% hist3 stuff

% $$$ lfp = Trial.lfp.copy;
% $$$ lfp.load(Trial,[69,80]);
% $$$ lfp.resample(xyz);
% $$$ wlfp = WhitenSignal(lfp.data,[],1);
% $$$ [ys,fs,ts] = mtcsdglong(wlfp(:,2),2^8,lfp.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
% $$$ 
% $$$ fncp = zeros(xyz.size(1),size(ys,2));
% $$$ tncp = zeros(xyz.size(1),1);
% $$$ fncp((2^6+1):(size(ys)+2^6),:) = ys;
% $$$ tncp((2^6+1):(size(ts)+2^6),:) = ts;
% $$$ fncp = MTADlfp('data',fncp,'sampleRate',xyz.sampleRate);
% $$$ %fncp.data = (fncp.data./repmat(sum(fncp.data),size(fncp.data,1),1));
% $$$ 
% $$$ tpow = MTADlfp('data',log10(mean(fncp(:,fs>=6&fs<=12),2)./mean(fncp(:,fs<4),2)),'sampleRate',xyz.sampleRate);;
% $$$ %tpow.filter(gtwin(1.2,tpow.sampleRate));
% $$$ 
% $$$ % $$$ figure,
% $$$ % $$$ sp(1) = subplot(211);
% $$$ % $$$ plot(tncp(nniz(tncp)),tpow(nniz(tncp)))
% $$$ % $$$ sp(2) = subplot(212);
% $$$ % $$$ imagesc(tncp(nniz(tncp)),fs,log10(fncp.data(nniz(tncp),:)')),axis xy
% $$$ % $$$ linkaxes(sp,'x');
% $$$ 
% $$$ 
% $$$ spowa = MTADlfp('data',nanmean(rhm(:,fs>6&fs<12),2),'sampleRate',Trial.xyz.sampleRate);
% $$$ spd = xyz.vel(marker,[1,2]);
% $$$ spd.data = log10(spd.data);
% $$$ 
% $$$ 
% $$$ sind = true(xyz.size(1),1);
% $$$ sind = Trial.stc{'w'};
% $$$ 
% $$$ %varA1 = log10(xyz(sind,7,3));
% $$$ varA1 = spd(sind);
% $$$ varA2 = log10(xyz(sind,7,3));
% $$$ %varA2 = spowa(sind);
% $$$ varA3 = spowa(sind);
% $$$ %varA3 = tpow(sind);
% $$$ 
% $$$ vind = nniz(varA1)&nniz(varA2)&nniz(varA3);
% $$$ 
% $$$ varA1 = varA1(vind);
% $$$ varA2 = varA2(vind);
% $$$ varA3 = varA3(vind);
% $$$ 
% $$$ vA1l = prctile(varA1,[2,98]);
% $$$ vA2l = prctile(varA2,[2,98]);
% $$$ vA3l = prctile(varA3,[2,98]);
% $$$ 
% $$$ vA1e = linspace(vA1l(1),vA1l(2),30);
% $$$ vA2e = linspace(vA2l(1),vA2l(2),30);
% $$$ vA3e = linspace(vA3l(1),vA3l(2),30);
% $$$ 
% $$$ vA12N = hist2([varA1,varA2],vA1e,vA2e);
% $$$ 
% $$$ [vA1c, vA1i] = histc(varA1,vA1e);
% $$$ [vA2c, vA2i] = histc(varA2,vA2e);
% $$$ 
% $$$ vA12i = [vA1i,vA2i];
% $$$ vA12in = vA12i(nniz(vA12i),:);
% $$$ varA3n = varA3(nniz(vA12i));
% $$$ 
% $$$ A = nan(29,29);
% $$$ A = accumarray(vA12in,varA3n,[29,29],@mean,nan);
% $$$ B = accumarray(vA12in,varA3n,[29,29],@std,nan);
% $$$ A(vA12N<30)=nan;
% $$$ B(vA12N<30)=nan;
% $$$ 
% $$$ figH = figure(3757),
% $$$ subplot(131);
% $$$ imagescnan({vA1e,vA2e,A'},prctile(A(nniz(A(:))),[5,95]),false,true,[0,0,0]),axis xy,
% $$$ title(['Mean NCP PSD']);
% $$$ xlabel('Log10 Body Speed (cm/s)');
% $$$ ylabel('Log10 Head Height (mm)');
% $$$ 
% $$$ subplot(132)
% $$$ imagescnan({vA1e,vA2e,B'},prctile(B(nniz(B(:))),[5,95]),false,true,[0,0,0]),axis xy,
% $$$ title(['Std NCP PSD']);
% $$$ xlabel('Log10 Body Speed (cm/s)');
% $$$ ylabel('Log10 Head Height (mm)');
% $$$ 
% $$$ subplot(133)
% $$$ ABN = A'./B';
% $$$ imagescnan({vA1e,vA2e,ABN},prctile(ABN(nniz(ABN(:))),[5,95]),false,true,[0,0,0]),axis xy,
% $$$ title(['SNR NCP PSD']);
% $$$ xlabel('Log10 Body Speed (cm/s)');
% $$$ ylabel('Log10 Head Height (mm)');
% $$$ 
% $$$ suptitle([Trial.filebase ' : ' Trial.stc{s}.label])
% $$$ set(figH,'pos',[14,325,1581,420]);
% $$$ 
% $$$ reportfig(fullfile(Trial.path.project,'figures'),figH, ...
% $$$           'meanNCP_on_jpdf_bv_zh',[],[Trial.filebase ' :' Trial.stc{s}.label],200)

% $$$ %figure,
% $$$ %subplot(133),imagescnan({vedges,sedges,A'},[-5,-3],[],1,[0,0,0]),axis xy,
% $$$ %clf
% $$$ %imagescnan({vedges,sedges,A'},[1.2,1.9],[],1,[0,0,0]),axis xy,
% $$$ %subplot(122),imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ % CA1 LM 81;
% $$$ % DG  G  85;
% $$$ % CA3 ?  95;
% $$$ chan = find(chans == 71);
% $$$ numIter = 10000;
% $$$ %tpow = log10(mean(yld(aind,fl>6&fl<12,chan),2));
% $$$ %tpow = log10(mean(yld(aind,fl<=4,chan),2));
% $$$ tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
% $$$ %tpow = log10(mean(yld(aind,fh>50&fh<80,chan),2));
% $$$ 
% $$$ B=nan(vbins,sbins,numIter);
% $$$ A=nan(vbins,sbins,numIter);
% $$$ %S=nan(vbins,sbins,numIter);
% $$$ tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$ B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
% $$$ A(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
% $$$ %S(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
% $$$ for i = 2:numIter,
% $$$ tpow = tpow(randperm(numel(tpow)));
% $$$ tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$ A(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
% $$$ %S(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
% $$$ end
% $$$ 
% $$$ AS = sort(A,3);
% $$$ P = 1./sum(repmat(A(:,:,1),[1,1,numIter])>A,3);
% $$$ P(isinf(P)) = nan;
% $$$ 
% $$$ SIG = P<=0.0002;
% $$$ ASIG = A; 
% $$$ ASIG(~SIG)=nan;
% $$$ ASIG(B<10)=nan;
% $$$ Aclims = [prctile(ASIG(~isnan(ASIG)),5),prctile(ASIG(~isnan(ASIG)),95)];
% $$$ 
% $$$ figure
% $$$ 
% $$$ subplot(131),imagescnan({vedges,sppedges,B'./yad.sampleRate},[],[],1,[0,0,0]),axis xy,
% $$$ title('Occupancy in seconds')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ subplot(132),imagescnan({vedges,sedges,   A(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
% $$$ title('Mean Power 1-4Hz given 10Hz Osc. Power dist(SU,HB) VS Vel(HF)')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ subplot(133),imagescnan({vedges,sedges,ASIG(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
% $$$ title('P<0.0002 and Occupancy > 1.33 seconds')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ flim=[1,4;6,12;20,27;30,40];
% $$$ %flim=[40,60;60,80;80,100;100,120];
% $$$ %mychans = [71,73,81,85,95];
% $$$ mychans = 1:2:8;
% $$$ B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
% $$$ figure
% $$$ for c = 1:numel(mychans),
% $$$     for i = flim',
% $$$         subplot2(numel(mychans),size(flim,1),c,find(i(1)==flim(:,1)));
% $$$         tpow = log10(mean(yld(aind,fl>i(1)&fl<i(2),find(chans==mychans(c))),2));
% $$$         %tpow = log10(mean(yld(aind,fh>i(1)&fh<i(2),find(chans==mychans(c))),2));
% $$$         tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$         AFB = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
% $$$         AFB(B<5)=nan;
% $$$         AFBclims = [prctile(AFB(~isnan(AFB)),5),prctile(AFB(~isnan(AFB)),95)];
% $$$         imagescnan({vedges,sedges,AFB'},AFBclims,[],1,[0,0,0]);
% $$$         axis xy,
% $$$         ticks_lin2log(gca,'x')
% $$$         title(['C: ' num2str(mychans(c)) 'Mean P(' num2str(i(1)) '-' num2str(i(2)) ')'])
% $$$     end
% $$$ end
% $$$ 
% $$$ text(.1,.1,['Mean P(' num2str(flim(1)) '-' num2str(flim(end)) ') given 10Hz Osc. Power dist(SU,HB) VS Vel(SL)'])
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ 
% $$$ 
% $$$ A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@nanmean,nan);
% $$$ figure,imagescnan({vedges,sedges,clip(A,0,.015)'},[],[],1,[0,0,0]),axis xy,
% $$$ % $$$ 
% $$$ % $$$ %tbp_phase.resample(yad);
% $$$ % $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_median,nan);
% $$$ % $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_var,nan);
% $$$ % $$$ figure,imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ 
% $$$ % $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4|(fl>12&fl<18),chan),2);
% $$$ % $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4,chan),2);
