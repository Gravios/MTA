
Session = MTATrial.validate('jg05-20120317.cof.all');

filebase = fullfile(Session.spath,Session.name);
channel = 72;
freqMethod = 'spect';
phaseMethod = 'hilbert';
filterType = 'adapt'; %'butt','cheb','mt'
overwrite = true;
signalType = 'lfp';
freqRange = [5,13];

[thPh, thAmp, thFr] = ThetaParams(filebase, channel, overwrite, freqRange, phaseMethod, freqMethod, filterType, ...
                                  signalType);


channel = 86;
freqMethod = 'spect';
phaseMethod = 'hilbert';
filterType = 'adapt'; %'butt','cheb','mt'
overwrite = true;
signalType = 'lfp';
freqRange = [5,13];

[thPhLM, thAmpLM, thFrLM] = ThetaParams(filebase,  channel,     overwrite, ...
                                        freqRange, phaseMethod, freqMethod, filterType, ...
                                        signalType);

lfp = Session.load('lfp',72);







%lfp = Session.load('lfp',[72,86]);
lfp = LoadBinary([filebase,'.lfp'],[72,86],96,[],[],[])';

freq_range = [5,13];
tbp = lfp;
tbp = ButFilter(tbp,3,[5,13]./(1250./2),'bandpass');
tbp_hilbert = zeros(size(lfp));
tbp_hilbert(nniz(lfp),:) = hilbert(tbp(nniz(lfp),:));
tbp_phase = phase(tbp_hilbert);




phz = lfp.phase([5,13]);    


figure,
timeLfp = (1:size(thPh,1))./1250;
sp = gobjects([0,1]);
sp(end+1)=subplot(5,2,[1,2]);
hold('on');
plot(timeLfp,log10(thAmp));
plot(timeLfp,log10(thAmpLM));
title('log10(thetaPower)');
sp(end+1)=subplot(5,2,[3,4]);
hold('on');
plot(timeLfp,thPh);
plot(timeLfp,thPhLM);
title('Full Session - labbox:TF:ThetaParam: theta phase CA1pyr and CA1lm');
sp(end+1)=subplot(5,2,[5,6]);
hold('on');
plot(timeLfp,[tbp_phase(:,1)]);
plot(timeLfp,[tbp_phase(:,2)]);
title('Full Session - MTA:MTAData:phase: theta phase CA1pyr and CA1lm');
sp(end+1)=subplot(5,2,[7,8]);
hold('on');
plot(linspace([Session.sync.data([1,end]),size(phz,1)]),[phz(:,1)]);
plot(linspace([Session.sync.data([1,end]),size(phz,1)]),[phz(:,2)]);
title('Synced Trial - MTAData:phase: theta phase CA1pyr and CA1lm');
linkaxes(sp,'x');
ind = thAmp>100;
sp(end+1)=subplot(5,2,[9]);
rose(circ_dist(tbp_phase(ind,1),thPh(ind,1)),1000);
title({'phase difference between','labbox:TF:ThetaParam and MTA:MTAData:phase','CA1pyr'});
sp(end+1)=subplot(5,2,[10]);
rose(circ_dist(tbp_phase(ind,2),thPhLM(ind,1)),1000);
title({'phase difference between','labbox:TF:ThetaParam and MTA:MTAData:phase','CA1lm'});
suptitle('Comparison between MTAData:phase and labbox:TF:ThetaParams');
