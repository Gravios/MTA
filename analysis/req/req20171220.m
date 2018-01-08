

sessionList = get_session_list('MjgER2016');

t = 1;

Trial = MTATrial.validate(sessionList(t));
xyz = preproc_xyz(Trial);

units = select_placefields(Trial);
pft = pfs_2d_theta(Trial);

drz = compute_drz(Trial, units, pft);
ddz = compute_ddz(Trial, units, pft);

%lfp = Trial.load('lfp',[Trial.ephys.electrode(:).CA1pyrThetaRef]);
lfp = Trial.load('lfp',sessionList(t).thetaRef);

phz = lfp.phase([6,12]);


spk = Trial.spk.copy();
spk.create(Trial,xyz.sampleRate,'theta-groom-sit',units,'deburst');    
spk.create(Trial,xyz.sampleRate,'loc&theta',units,'deburst');    

unit = units(1);
[P,phzStats] = MjgER2016_phasePrecession(Trial,drz,phz,spk,unit);