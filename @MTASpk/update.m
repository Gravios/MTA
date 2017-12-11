function Spk = update(Spk,Session)
% function Spk = update(Spk,Session)
% 
% Reload and map units from clu files into MTASpk object
%
Session  = MTASession.validate(Session);
Session.spk.create(Session);
Session.save;
compute_neuron_quality(Session,'overwrite',true);
