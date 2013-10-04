
Session = MTATrial('jg05-20120317');
Session = Session.load_ang(0);

hb = Session.transformOrigin('head_back','head_front',{'head_left','head_right'});


wag = circ_dist(Session.ang(:,1,4,1),Session.ang(:,1,2,1));
wig = circ_dist(Session.ang(:,4,1,1),Session.ang(:,4,3,1));

wag(isnan(wag))=0;
wig(isnan(wig))=0;

wigwag = WhitenSignal(wag-wig);

rol = hb.roll;
rol(isnan(rol))=0;

wol = WhitenSignal(rol);

[y,f,t,phi,fStats] = mtchglong([wol,wigwag],2^9,Session.xyzSampleRate,2^8,2^8-1,[],[],[],[1,60]);