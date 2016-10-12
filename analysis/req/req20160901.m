

%% Explore the use of marker subesets

MTAstartup('vr_exp');
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

sessionList = 'Ed10VR';

%% Load and preprocesses data

S = get_session_list(sessionList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');

Trial = MTATrial.validate(S(1));




xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,50);
ang = create(MTADang,Trial,xyz);

figure,
plot(circ_dist(ang(:,1,2,2),ang(:,1,3,2)))
hold on
plot(nunity(xyz(:,1,3)))

figure,plot(nunity(ang(:,2,3,3)))



wang.data = circ_dist(ang(:,1,2,2),ang(:,1,3,2));
wang.filter('ButFilter',3,[1,8],'bandpass');
figure,plot(wang.data)

sang = xyz.copy;
sang.data = circ_dist(ang(:,1,4,1),ang(:,1,5,1));
sang.filter('ButFilter',3,[1,8],'bandpass');
figure,plot(sang.data)


figure,hold on
plot(wang.data)
plot(-sang.data)

dang = xyz.copy;
dang.data = minus(ang(:,4,2,3),ang(:,4,3,3));
dang.filter('ButFilter',3,[1,8],'bandpass');
figure
plot(dang.data)




figure,hold on
plot(nunity(ang(:,1,2,3))/10)
plot(circ_dist(ang(:,'spine_lower','spine_upper',1),...
               circshift(ang(:,'spine_lower','spine_upper',1),-10)))
