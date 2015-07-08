function out  = fet_x_head_spine(Session)


ang = create(Session.ang.copy,Session,Session.load('xyz').filter(gtwin(.01,Session.xyz.sampleRate)));

out = [ButFilter(circ_dist(ang(:,1,4,1),ang(:,2,4,1)),3,1/(ang.sampleRate/2),'high'),...
       circ_dist(ang(:,1,4,1),ang(:,2,4,1)),...
       circ_dist(ang(:,5,7,1),ang(:,1,4,1))./2.2];
       

