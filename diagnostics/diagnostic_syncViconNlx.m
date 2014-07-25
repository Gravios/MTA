
s = MTASession('Ed03-20140624','cof',true,'0x0002');
figure,plot(s.xyz(:,7,3))
s.xyz.size
s.ang.size
clear

s = MTASession('Ed03-20140624');
s.load('xyz');
s.load('ang');
s.xyz.size
s.ang.size
hold on,plot(s.xyz(:,3,3),'r')

