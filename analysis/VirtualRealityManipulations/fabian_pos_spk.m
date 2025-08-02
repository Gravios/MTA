%% relinking data ------------------------------------------------------
dPaths.xyz   = fullfile('/storage/gravio/data/processed/xyz/FS04');   %|
dPaths.ephys = fullfile('/storage/gravio/data/processed/ephys/FS04'); %|
link_session('FS04-20210322a',dPaths);                                %|
%-----------------------------------------------------------------------


Trial = MTATrial.validate('FS04-20210322a.vrr.all');
rat   = Trial.load('subject','FS04');
rat   = Trial.load('subject','FS04_AC');
Arena = Trial.load('subject','Arena');


spk = Trial.load('spk',rat.sampleRate,'theta',[],'');

state = [Trial.stc{'theta'}];



hfig = figure();
for u = spk.map(1:end,1)'
    clf(hfig);
    hold('on');
    res = spk(u);
    plot(rat(state,'Head',1),rat(state,'Head',2),'.');
    plot(rat(res,'Head',1),rat(res,'Head',2),'.r');
    waitforbuttonpress();
end