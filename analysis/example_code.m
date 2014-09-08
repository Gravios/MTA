%% How to 

Session = 'Ed10-20140817';
overwrite = true;

xyz_path = '/gpfs01/sirota/homes/eduardo/data/xyz';
nlx_path = '/gpfs01/sirota/homes/eduardo/data/rawnlx';
linkSession(Session,xyz_path,nlx_path);



Session = MTASession('Ed10-20140817','cof',true,'0x0002');

plot(Session.xyz(:,'head_front',3));

Trial = QuickTrialSetup(Session,'dark',[],[4:8]);
%Trial = QuickTrialSetup(Session);

clear


Trial = MTATrial('Ed10-20140817');
xyz = Trial.xyz.copy

% figure,plot(xyz(Trial.stc{'r'},'head_front',3))

lfp_ncp = Trial.lfp.copy;
lfp_ncp.load(Trial,65);

lfp_olf = Trial.lfp.copy;
lfp_olf.load(Trial,33:64);


