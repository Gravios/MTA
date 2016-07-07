function Trial = labelAuxBhv(Trial,varargin)
[Stc,angThresh,overwrite] = DefaultArgs(varargin,{[],-.45,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
end


sempty = isempty(Stc{'h'})||isempty(Stc{'l'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('h')) = [];
        Stc.states(Stc.gsi('l')) = [];
    end
    
    fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
    if ~exist(fpath,'file'),
        bhv_rhm_distrb(Trial);
    end

    afig = hgload(fpath);
    rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','bhv_rhm_distrb-hangle'));
    rhm_distrb = rhm_distrb(1);
    figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
    mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
    mrhmp(mrhmp==0) = nan;
    rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
        +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
    angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
    delete(afig);
    % Split walking state into low and high walk
    wind = Stc{'w'}.copy;
    lang = Stc{'w'}.copy;
    hang = Stc{'w'}.copy;    

    ang = create(MTADang,Trial,Trial.load('xyz'));

    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,round(.1/ang.sampleRate));
    lang.label = 'lwalk';
    lang.key = 'l';
    lang.filename = [Trial.filebase '.sst.lwalk.l.mat'];
    Stc.states{end+1} = lang&wind.data;

    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(.1/ang.sampleRate));
    hang.label = 'hwalk';
    hang.key = 'h';
    hang.filename = [Trial.filebase '.sst.hwalk.h.mat'];
    Stc.states{end+1} = hang&wind.data;
    
end

Stc.save(1);
Trial.stc = Stc;
Trial.save;


% $$$ sempty = isempty(Stc{'h'})||isempty(Stc{'l'});
% $$$ if sempty||overwrite,
% $$$     [rhm,fs,ts] = fet_rhm(Trial,[],'mtcsdglong',true);
% $$$     xyz = Trial.xyz.copy;xyz.load(Trial);xyz.filter('ButFilter',3,4);
% $$$     ang = create(MTADang,Trial,xyz);
% $$$ 
% $$$     nrhm = rhm.copy;
% $$$     nrhm.data = log10(nrhm.data);
% $$$     nrhm.data(nrhm.data<-9) = nan;
% $$$ 
% $$$     rhmpow =median(nrhm(:,fs>5&fs<14),2);
% $$$     rhmpow = MTADlfp('data',rhmpow,'sampleRate',nrhm.sampleRate);
% $$$     ang.resample(rhmpow);
% $$$     xyz.resample(rhmpow);
% $$$ 
% $$$     rhmlims = linspace(-6,-2,40);
% $$$     anglims = linspace(-1.4,1.4,40);
% $$$ 
% $$$     bper = Stc{'w',rhmpow.sampleRate}.copy;
% $$$     bper.cast('TimeSeries');
% $$$     if bper.size(1)-ang.size(1)<0,bper.data=[bper.data;false([abs(bper.size(1)-ang.size(1)),1])];end
% $$$     xaind = nniz(ang(:,5,7,2))&nniz(rhmpow.data)&bper(1:ang.size(1));
% $$$     fet = [ang(xaind,5,7,2),rhmpow(xaind)];
% $$$     [bsts,bhmm,bdcd] = gausshmm(fet,2);
% $$$ 
% $$$     fasts = zeros([ang.size(1),1]);
% $$$     fasts(xaind) = bsts;
% $$$ 
% $$$ 
% $$$ 
% $$$     if mean(ang(fasts==1,5,7,2))<mean(ang(fasts==2,5,7,2))
% $$$         g = 2;l = 1;
% $$$     else
% $$$         g = 1;l = 2;
% $$$     end
% $$$     Stc.states(Stc.gsi('h')) = [];
% $$$     Stc.addState(Trial.spath,...
% $$$                  Trial.filebase,...
% $$$                  ThreshCross(fasts==g,.5,1),...                   
% $$$                  rhmpow.sampleRate,...
% $$$                  Trial.sync.copy,...
% $$$                  Trial.sync.data(1),...
% $$$                  'hswalk','h','TimePeriods');
% $$$     Stc.states{Stc.gsi('h')} = Stc{'h',Trial.xyz.sampleRate}&Stc{'w'};
% $$$     Stc.states{Stc.gsi('h')} = Stc{'h'}+[1/Trial.xyz.sampleRate,-1/Trial.xyz.sampleRate];
% $$$     %Stc.states{Stc.gsi('h')}.cast('TimePeriods');
% $$$ 
% $$$     Stc.states(Stc.gsi('l')) = [];
% $$$     Stc.addState(Trial.spath,...
% $$$                  Trial.filebase,...
% $$$                  ThreshCross(fasts==l,.5,1),...
% $$$                  rhmpow.sampleRate,...
% $$$                  Trial.sync.copy,...
% $$$                  Trial.sync.data(1),...
% $$$                  'lswalk','l','TimePeriods');
% $$$     %Stc.states{Stc.gsi('l')}.cast('TimePeriods');
% $$$     Stc.states{Stc.gsi('l')} = Stc{'l',Trial.xyz.sampleRate}&Stc{'w'};
% $$$ end






