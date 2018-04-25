
Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial = MTATrial.validate('jg04-20120213.cof.all');
Trial = MTATrial.validate('jg04-20120210.cof.all');

s = MTASession.validate(Trial.filebase);
s.spk.create(s);
s.save();
Trial = MTATrial.validate(Trial.filebase);
%compute_neuron_quality(Trial,[],[],[],true);


Trial.load('nq');
nq = Trial.nq;

units = select_units(Trial);


units = select_placefields(Trial,[],false);

pft = pfs_2d_theta(Trial,[],false,true,1);
%pft.purge_savefile();

figure,
for unit = pft.data.clu(:)',
%for unit = units,
    plot(pft,unit,1,true,[],true);
    title(num2str(unit))
    waitforbuttonpress();
end



[dfs,dfsTags] = req20180123_ver5(Trial,[],'7',true);