function feature = fet_eduardo(Session)

xyz = Session.load('xyz');

%feature = xyz.copy;
%feature.label = 'error';
%feature.key   = 'e';

feature = xyz(:,'head_front',3);