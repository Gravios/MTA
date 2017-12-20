


Trial = MTATrial.validate('er01-20110719.cof.all');

units = select_placefields(Trial);
pft = pfs_2d_theta(Trial);

drz = compute_drz(Trial,units,pft);



[P,phzStats] = MjgER2016_phasePrecession(Trial,drz,phz,spk,unit,state);