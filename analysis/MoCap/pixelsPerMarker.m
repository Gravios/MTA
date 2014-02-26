
%% Mcam
hpx = 953;
vpx = 881;
msize = 8 %mm
fov = 41;

%%18w
hpx = 1664;
vpx = 1088;
msize = 8 %mm
fov = 70;


%%p18w
hpx = 1664;
vpx = 1088;
msize = 4 %mm
fov = 67;


%%p41
hpx=2048;
vpx=2048;
msize=8;
fov = 51;


fov = deg2rad(fov);
r = [.4:.01:4]';

c = sqrt((r-r.*cos(fov)).^2+(r.*sin(fov)).^2);

hppm = msize./(c./hpx.*1000);
vppm = msize./(c./vpx.*1000);


mazer = [.5:.1:4]';

dr = (2.*r./sin((pi-fov)./2)).*sqrt(1/((1-cos(fov)).^2+sin(fov).^2)).*cos(fov./2);


figure
plot([r,r,r,r],[hppm,vppm,c,dr])



