% req20170926
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: examination of new RectFilter and usefulness in fet_rhm
%  Bugs: NA


Trial = MTATrial.validate('jg05-20120317.cof.all');

xyz = Trial.load('xyz');
fxyz = xyz.copy();
fxyz.filter('RectFilter');
ang = create(MTADang,Trial,xyz);fang = create(MTADang,Trial,fxyz);

plot_features_with_stc(MTADfet.encapsulate(Trial,[circ_dist(circshift(fang(:,2,4,1),-1),circshift(fang(:,2,4,1),1)),circ_dist(circshift(fang(:,5,7,1),-1),circshift(fang(:,5,7,1),1))],xyz.sampleRate,'azvel','azvel','v'),Trial.load('stc','hand_labeled_rev3_jg'))


plot_features_with_stc(MTADfet.encapsulate(Trial,[circ_dist(circshift(fang(:,2,4,2),-1),circshift(fang(:,2,4,2),1)),circ_dist(circshift(fang(:,5,7,2),-1),circshift(fang(:,5,7,2),1))],xyz.sampleRate,'azvel','azvel','v'),Trial.load('stc','hand_labeled_rev3_jg'))


hbfet = fet_HB_pitchvel(Trial);
plot_features_with_stc(hbfet,Trial.load('stc','hand_labeled_rev3_jg'));

hbfet = fet_HB_angvel(Trial);
plot_features_with_stc(hbfet,Trial.load('stc','hand_labeled_rev3_jg'));


