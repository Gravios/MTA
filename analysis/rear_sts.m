function rear_sts(Session)

%Session = MTASession(Session);
%Session = MTASession('jg05-20120310',{'ang'});
Session = MTASession(Session,{'ang'});

rper = Session.Bhv.getState('rear').state;
rdur = diff(rper,1,2);
rfet = Session.ang(:,Session.Model.gmi('spine_middle'),Session.Model.gmi('spine_upper'),2).*Session.xyz(:,Session.Model.gmi('head_back'),3);
rfet(isnan(rfet))=-1;

rseg = GetSegs(rfet,round(rper-5*Session.xyzSampleRate),round(10*Session.xyzSampleRate),[]);
rseg = reshape(rseg,[],size(rper,1),2);

gri = find(max(rseg(1:round(size(rseg,1)/3),:,1))<20&rdur'>1.5*Session.xyzSampleRate);
griu = find(max(rseg(round(size(rseg,1).*.66):end,:,2))<20&rdur'>1.5*Session.xyzSampleRate);

tvec = [-1250,-500;-750,0;-375,375;0,750;500,1250;-1250,-500;-750,0;-375,375;0,750;500,1250];

rsts = {'pre_rear_onset','pre_peri_rear_onset' 'dia_rear_onset','post_peri_rear_onset', 'post_rear_onset';'pre_rear_offset','pre_peri_rear_offset','dia_rear_offset','post_peri_rear_offset','post_rear_offset'};
rtmp = round(reshape(tvec,5,2,2)./1000.*Session.lfpSampleRate);
for b = 1:2,
    for state = 1:length(rsts),
        trsts = round((repmat(sq(rtmp(state,b,:))',length(gri),1) + repmat(rper(gri,b),1,2)-1)./Session.xyzSampleRate.*Session.lfpSampleRate+Session.syncPeriods(1,1));
        msave([Session.spath.nlx Session.name '.sts.' rsts{b,state}],trsts);
    end
end



