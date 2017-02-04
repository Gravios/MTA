



durThresh = 1;
Bccg{s} = gen_bhv_ccg(Trial,'walk',durThresh);

figure,
for i = 1:100,
    subplot(211); cla;
    Bccg.plot(i,1);
    subplot(212); cla;
    Bccg.plot(i,2);
    pause(.2) 
end