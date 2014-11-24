function Trial = labelAllBhv(Trial)

Trial.stc.states = {};
try,Trial = labelBhv(Trial);end
try,Trial = labelAuxBhv(Trial,'overwrite',true);end
try,Trial.stc.states{end+1} = theta(Trial);end
Trial.stc.save(1);
Trial.save;