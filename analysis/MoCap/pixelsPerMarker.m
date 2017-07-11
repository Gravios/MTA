
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


%%p17w
cam(3).hpx = 1664;
cam(3).vpx = 1088;
cam(3).msize = 8; %mm
cam(3).fov = 70;


%%p41
cam(4).name = 'Prime 41';
cam(4).hpx=2048;
cam(4).vpx=2048;
cam(4).msize=5;
cam(4).fov = 51;

cam(5).name = 'Prime 13';
cam(5).hpx=1280;
cam(5).vpx=1024;
cam(5).msize= 5;
cam(5).fov = 56;



%for i=1:numel(c)
%end

i = 4;

cam(i).fov = deg2rad(cam(i).fov);
r = [.4:.01:4]';
c = sqrt((r-r.*cos(cam(i).fov)).^2+(r.*sin(cam(i).fov)).^2);

hppm = cam(i).msize./(c./cam(i).hpx.*1000);
vppm = cam(i).msize./(c./cam(i).vpx.*1000);


mazer = [.5:.1:4]';

dr = (2.*r./sin((pi-cam(i).fov)./2)).*sqrt(1/((1-cos(cam(i).fov)).^2+sin(cam(i).fov).^2)).*cos(cam(i).fov./2);


figure
plot([r,r,r,r.*2],[hppm,vppm,c,dr])
legend(['Horizontal pixel count for the cross section of a ' num2str(cam(i).msize) ' mm marker given distance to camera'],...
       ['Vertical pixel count for the cross section of a ' num2str(cam(i).msize) ' mm marker given distance to camera'],...
       ['Maximum maze diameter given distance between maze center and camera (m)'],...
       ['Required distance between maze center and camera given maze diameter (m)'])
title({['Camera Name: ' cam(i).name],...
       ['Camera Res: ' num2str(cam(i).hpx) 'x' num2str(cam(i).vpx)  'pixels'],...
       ['Camera fov: ' num2str(rad2deg(cam(i).fov))]})
