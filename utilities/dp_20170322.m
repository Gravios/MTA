
sessionList = get_session_list('BHV_S4H5'); 
s = 6;%19
listFiles(sessionList(s).sessionName,'.trl.')
listFiles(sessionList(s).sessionName,'stc.NN0317')
listFiles(sessionList(s).sessionName,'.trb.')
listFiles(sessionList(s).sessionName,'3dssh')
PlotSessionErrors(MTATrial.validate(sessionList(s)));
labelBhv_NN(MTATrial.validate(sessionList(s)));

arrayfun(@(s) listFiles(s.sessionName,'stc.NN0317.'),sessionList,'UniformOutput',false)
