linkSession('RS0615-20151201S1',...
            '/storage/gravio/data/processed/xyz/RS0615/',...
            '/storage/ricardo/data/Freely_moving/RS0615/Neuralynx/');

s = MTASession('RS0615-20151201S1','sof',true,'0x0008','ephySystem','bob');

%Data = wrapBS(s,rec.corrsub.);
Data.save


chx = s.load('fet','cbs');

xyz = s.load('xyz');
chx.resample(xyz);

vxy = xyz.vel;



fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
ang = create(MTADang,s,xyz);


[rhm,fs,ts] = fet_rhm(s,[],'mtchglong',false);

figure,imagesc(ts,fs,log10(rhm.data')),caxis([-9,-3]),colormap jet, ...
    axis xy


fet = s.fet.encapsulate

name = 'head_back_height'; label = 'hb_z'; key = 'h';
fet = MTADfet.encapsulate(Trial,...
                         xyz(:,1,3),...
                         xyz.sampleRate,...
                         name,label,key);

[ys,fs,ts] = fet_spec(Trial,fet);

figure,imagesc(ts,fs,log10(ys.data')),caxis([-9,-3]),colormap jet, ...
    axis xy
