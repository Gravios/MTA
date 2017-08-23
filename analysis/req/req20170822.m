




% Removed other states where shake occured. 
% Final bug in label_bhv_shake eliminated.
Stc = cf(@(t) t.load('stc','msnn_ppsvd'),Trials);
for t = 1:numel(Trials),
aper = Stc{t}{'k'};
kper = Stc{t}{'k'};
cf(@(s,p) cf(@(sts,per) set(sts,'data',get(sts-bsxfun(@plus,per,[-1,1]),'data')),...
             s.states,repmat({p},[1,numel(s.states)])),...
   Stc(t),repmat({aper.data},[1,1]));
Stc{t}.states{7} = kper;
end
cf(@(s) s.save(1),Stc);
t = 3;
fetb = fet_bref(Trials{t});
figure,plot(fetb.data(:,17:2:21)); plotSTC(Stc{t},[],[],{'walk','rear','turn','pause','groom','sit','shake'},'brgcymk');


sessionList = 'hand_labeled';
Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
states = {'walk','rear','turn','pause','groom','sit','shake'};
stc = 'msnn_ppsvd';

labelingStats = compute_inter_stc_stats(sessionList,stc,states,119.881035);





StcHL = cf(@(t) t.load('stc'),Trials);
t = 1;
fetb = fet_bref(Trials{t});
figure,sp = [];
sp(end+1) = subplot2(4,1,1,1);plot(fetb.data(:,16:2:20)); 
sp(end+1) = subplot2(4,1,2,1);plot(fetb.data(:,17:2:21)); 
sp(end+1) = subplot2(4,1,3,1);plotSTC(Stc{t},[],[],{'walk','rear','turn','pause','groom','sit','shake'},'brgcymk');
sp(end+1) = subplot2(4,1,4,1);plotSTC(StcHL{t},[],[],{'walk','rear','turn','pause','groom','sit','shake'},'brgcymk');

linkaxes(sp,'x');


