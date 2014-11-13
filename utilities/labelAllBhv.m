function Trial = labelAllBhv(Trial)

try,Trial = labelBhv(Trial);end
try,Trial = labelAuxBhv(Trial);end
try,Trial.stc.states{end+1} = theta(Trial);end
Trial.stc.save(1);
Trial.save;