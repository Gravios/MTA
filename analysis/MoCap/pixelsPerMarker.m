
%% Mcam
cam(1).hpx = 953;
cam(1).vpx = 881;
cam(1).msize = 8; %mm
cam(1).fov = 41;

%%18w
cam(2).hpx = 1664;
cam(2).vpx = 1088;
cam(2).msize = 8; %mm
cam(2).fov = 51;


%%p18w
cam(3).hpx = 1664;
cam(3).vpx = 1088;
cam(3).msize = 8; %mm
cam(3).fov = 70;


%%p41
cam(4).hpx=2048;
cam(4).vpx=2048;
cam(4).msize=8;
cam(4).fov = 51;




for i=1:numel(c)
end

i = 3

cam(i).fov = deg2rad(cam(i).fov);
r = [.4:.01:4]';
c = sqrt((r-r.*cos(cam(i).fov)).^2+(r.*sin(cam(i).fov)).^2);

hppm = cam(i).msize./(c./cam(i).hpx.*1000);
vppm = cam(i).msize./(c./cam(i).vpx.*1000);


mazer = [.5:.1:4]';

dr = (2.*r./sin((pi-cam(i).fov)./2)).*sqrt(1/((1-cos(cam(i).fov)).^2+sin(cam(i).fov).^2)).*cos(cam(i).fov./2);


figure
plot([r,r,r,r],[hppm,vppm,c,dr])


camera(1).pos = 
