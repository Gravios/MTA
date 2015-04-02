function out  = fet_p_head_spine(Session)


ang = create(Session.ang.copy,Session,Session.load('xyz').filter(gtwin(.01,Session.xyz.sampleRate)));

out = [ang(:,3,4,2),ang(:,5,7,2)];