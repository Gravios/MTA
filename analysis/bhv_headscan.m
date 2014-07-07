%bhv_headscan.m
%quick look at head scanning motions

Trial = MTATrial('jg05-20120317');
Trial.load('xyz');
Trial.xyz.filter(gausswin(31)./sum(gausswin(31)));
ang = Trial.ang.copy;
ang.create(Trial,Trial.xyz);
figure,plot(circ_dist(ang(:,4,5,1),ang(:,5,7,1)));
hold on,plot(diff(Filter0(gausswin(31)./sum(gausswin(31)),circ_dist(ang(:,4,5,1),ang(:,5,7,1)))).*10,'r');
title('Neck Joint angle and angular speed')
xlabel time