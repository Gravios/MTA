function units = select_placefields(Trial)

Trial.load('nq');
pft = pfs_2d_theta(Trial); 
mrt = pft.maxRate;
units = select_units(Trial,18,'pyr');
%units = units(Trial.nq.SNR(units)>1);
units = units(mrt(pft.data.clu(units))>0.8);
units = units(pft.data.si(units)>0.5);
units = units(pft.data.spar(units)<0.5);

