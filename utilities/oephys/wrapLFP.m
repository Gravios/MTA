function Data = wrapLFP(Session,data)
Data = MTADfet(Session.spath,...
               Session.filebase,...
               data,...
               1250,...
               Session.sync.copy,...
               Session.sync.data(1),...
               [],...
               'TimeSeries',...
               [],...
               'Local Field Potential',...
               'lfp',...
               'l');

