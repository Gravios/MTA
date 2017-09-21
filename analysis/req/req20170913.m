% req20170811 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Single Step Analysis
%  Bugs: NA

Trial = MTATrial.validate('Ed05-20140529.ont.all');

stc = Trial.load('stc','hand_labeled_rev1_Ed');

features = fet_bref(Trial);
ffet = features.copy();
ffet.filter('ButFilter',3,[2.5],'low');


xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,[2.5],'low');
vxy = xyz.vel(1,[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

figure,
subplot2(3,1,[1,2],1);plot(features(:,17));
subplot2(3,1,3,1);plotSTC(stc);
linkaxes(findobj(gcf,'Type','axes'),'x');


[steps,sval] = LocalMinima(-abs(features(:,17)),10,-4);

stdist = cf(@(st,s)  sum(WithinRanges(st,s.data)), repmat({steps},1,numel(stc.states)),stc.states);


fphase = features.phase([1.2,6]);

features.filter('ButFilter',3,[1,6],'bandpass');
hphase = features.phase([1,6]);

fh = hilbert(features.data(:,[17,25]));
fppc = zeros(size(steps));
sppc = zeros(size(steps));
hppc = zeros(size(steps));
hppa = zeros(size(steps));
for i = 1:numel(steps),
    %fppc(i) = PPC(fphase(steps(i),[17,19,21,23]));
    fppc(i) = circ_dist(fphase(steps(i),[17]),fphase(steps(i),[19]));
    sppc(i) = circ_dist(fphase(steps(i),[17]),fphase(steps(i),[21]));
    hppc(i) = circ_dist(fphase(steps(i),[17]),fphase(steps(i),[25]));
    hppa(i) = abs(fh(steps(i),2));
end



steps = LocalMinima(-abs(fh(:,2)),10,-4);
fppc = circ_dist(fphase(steps,[17]),fphase(steps,[19]));
sppc = circ_dist(fphase(steps,[17]),fphase(steps,[21]));
hppc = circ_dist(fphase(steps,[17]),fphase(steps,[25]));
hppad = circ_dist(angle(fh(steps,[1])),angle(fh(steps,[2])));

figure,
plot(abs(fh(:,2)))

figure,rose(hppc,100)
figure,rose(sppc,100)
figure,rose(fppc,100)

figure,
subplot(131);
plot(fppc,log10(abs(fh(steps,2))),'.')
subplot(132);
plot(sppc,log10(abs(fh(steps,2))),'.')
subplot(133);
plot(hppc,log10(abs(fh(steps,2))),'.')

figure,
plot(angle(fh(steps,2)),log10(abs(fh(steps,2))),'.')

figure,
plot(hppad,log10(abs(fh(steps,2))),'.')




figure,plot(fppc,log10(hppa),'.')
figure,plot(hppc,log10(hppa),'.')


figure(); hold('on');
subplot(221);plot(sppc,fppc,'.');grid('on')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(sppc(ind),fppc(ind),'.r');grid('on')
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(sppc(ind),fppc(ind),'.c');grid('on')
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(sppc(ind),fppc(ind),'.g');grid('on')
%ForAllSubplots('xlim([0.5,1.5]),ylim([-0.2,1])')
ForAllSubplots('xlim([-pi,pi]),ylim([-pi,pi])')

figure(); hold('on');
subplot(221);plot(log10(-sval),fppc,'.')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(log10(-sval(ind)),fppc(ind),'.r');
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(log10(-sval(ind)),fppc(ind),'.c');
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(log10(-sval(ind)),fppc(ind),'.g');
%ForAllSubplots('xlim([0.5,1.5]),ylim([-0.2,1])')
ForAllSubplots('xlim([0.5,1.5]),ylim([-2,2])')

figure(); hold('on');
subplot(221);plot(vxy(steps),fppc,'.');grid('on')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(vxy(steps(ind)),fppc(ind),'.r');grid('on')
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(vxy(steps(ind)),fppc(ind),'.c');grid('on')
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(vxy(steps(ind)),fppc(ind),'.g');grid('on')
%ForAllSubplots('xlim([0.5,1.5]),ylim([-0.2,1])')
ForAllSubplots('xlim([-pi,pi]),ylim([-pi,pi])')

figure(); hold('on');
subplot(221);plot(vxy(steps),hppc,'.');grid('on')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(vxy(steps(ind)),hppc(ind),'.r');grid('on')
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(vxy(steps(ind)),hppc(ind),'.c');grid('on')
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(vxy(steps(ind)),hppc(ind),'.g');grid('on')
%ForAllSubplots('xlim([0.5,1.5]),ylim([-0.2,1])')
ForAllSubplots('xlim([-pi,pi]),ylim([-pi,pi])')



figure(); hold('on');
subplot(221);plot(features(steps,16),log10(-sval),'.');grid('on')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(features(steps(ind),16),log10(-sval(ind)),'.r');grid('on')
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(features(steps(ind),16),log10(-sval(ind)),'.c');grid('on')
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(features(steps(ind),16),log10(-sval(ind)),'.g');grid('on')
linkaxes(findobj(gcf,'Type','axes'),'xy');


figure(); hold('on');
subplot(221);plot(ffet(steps,16),log10(-sval),'.');grid('on')
ind = WithinRanges(steps,stc{'w'}.data);
subplot(222);plot(ffet(steps(ind),16),log10(-sval(ind)),'.r');grid('on')
ind = WithinRanges(steps,stc{'p'}.data);
subplot(223);plot(ffet(steps(ind),16),log10(-sval(ind)),'.c');grid('on')
ind = WithinRanges(steps,stc{'n'}.data);
subplot(224);plot(ffet(steps(ind),16),log10(-sval(ind)),'.g');grid('on')
linkaxes(findobj(gcf,'Type','axes'),'xy');


sper = zeros([size(features,1),1]);
sper(steps) = 1;
sper = conv(sper,gausswin(120)./sum(gausswin(120)),'same');


figure,
subplot2(3,1,[1,2],1);plot(sper);
subplot2(3,1,3,1);plotSTC(stc);
linkaxes(findobj(gcf,'Type','axes'),'x');

