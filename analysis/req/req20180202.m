

Trial = MTATrial.validate('jg05-20120317.cof.all');
stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = preproc_xyz(Trial,'trb');
ang = create(MTADang,Trial,xyz);

pitchReferenceTrial = 'Ed05-20140529.ont.all';



pch = fet_HB_pitch(Trial);
map_to_reference_session(pch,Trial,pitchReferenceTrial);    


edx = linspace([-pi/2,pi/2,60]);
edy = edx;


figure,
ind = [stc{'x+p+r'}];
subplot(2,2,1);
hist2([pch(ind,1),pch(ind,3)],edx,edy);
xlabel('pitch SM_SU');
ylabel('pitch SU_HC');
subplot(2,2,2);
hist2([pch(ind,1),circ_dist(pch(ind,3),pch(ind,1))],edx,edy);
xlabel('pitch SM_SU');
ylabel('pitch diff(SM_SU,SU_HC)');
subplot(2,2,3);
hist2([pch(ind,1),pch(ind,3)],edx,edy);
subplot(2,2,4);
hist2([pch(ind,1),circ_dist(pch(ind,3),pch(ind,1))],edx,edy);
hold('on')




for s = [stc{'lloc'}.data]',
    ind = s(1):10:s(2);
    plot(pch(ind,1),circ_dist(pch(ind,3),pch(ind,1)),'.g','MarkerSize',1)
end
for s = [stc{'hloc'}.data]',
    ind = s(1):10:s(2);
    plot(pch(ind,1),circ_dist(pch(ind,3),pch(ind,1)),'.r','MarkerSize',1)
end




for s = stc{'r'}.data',
    ind = (s(1)-10):(s(1)+60);
    plot(pch(ind,1),circ_dist(pch(ind,3),pch(ind,1)),'.m','MarkerSize',1)
end
for s = stc{'r'}.data',
    ind = (s(2)-60):(s(2)+10);
    plot(pch(ind,1),circ_dist(pch(ind,3),pch(ind,1)),'.c','MarkerSize',1)
end

