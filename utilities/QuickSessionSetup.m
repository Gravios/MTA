function QuickSessionSetup(SessionParm)

if ischar(SessionParm), 
Sessions = SessionList(SessionParm);
else, Sessions = {SessionParm};end

for s = 1:numel(Sessions)
MTAstartup([],Sessions{s}{4});
Session = MTASession(Sessions{s}{1},Sessions{s}{2},true,Sessions{s}{5},Sessions{s}{6});
end