function rand_surrogate_segments(Session,state,nSegs,segLen)

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

stsp = Session.statePeriods('nrhp');

surrogateIndex = [];
for i = 1:size(stsp,1),
    surrogateIndex = cat(1,surrogateIndex,[stsp(i,1):stsp(i,2)]');
end

rpi = randi(size(surrogateIndex,1),nSegs,1);


save([Session.spath.nlx Session.name '.sts.rand_sur_seg