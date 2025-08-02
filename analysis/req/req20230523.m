% example data


MjgER2016_load_data();

Trial = Trials{20};

units = select_units(Trial,'pyr');

xyz = preproc_xyz(Trial,'trb');

spk = Trial.load('spk',xyz.sampleRate,'',units);

[states_matrix,states_labels] = stc2mat(Trial.stc,xyz,{'theta','walk','turn','pause','rear','groom','sit','loc','lloc','hloc','lpause','hpause'});
states_matrix = logical(states_matrix);

spikes_timestamps = spk.res;
spikes_clusterIds = spk.clu;
spikes_ids = unique(spikes_clusterIds);

position = sq(xyz(:,'hcom',[1,2]));
sample_rate = xyz.sampleRate;


save('/storage/gravio/data/example_data.mat','sample_rate','position','spikes_timestamps','spikes_clusterIds','spikes_ids','states_matrix','states_labels');


index = states_matrix(:,1) & ~any(states_matrix(:,[6,7]),2);
figure,
plot(position(index,1),position(index,2),'.')
hold('on');
mspk = spikes_timestamps(spikes_clusterIds==20);
spkIndex = mspk(index(mspk));
plot(position(spkIndex ,1),position(spkIndex,2),'.r');