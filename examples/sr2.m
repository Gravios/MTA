


MjgER2016_load_data();

Trial = Trials{20};
unitSubset = units{20};


Trial = MTATrial.validate('jg05-20120312.cof.all');




xyz = Trial.load('xyz','trb');

figure,
plot(xyz(:,'head_back',1),xyz(:,'head_back',2),'.');


fxyz = filter(copy(xyz),'ButFilter',3,2.5,'low');
vxy = vel(fxyz,{'spine_lower','head_front'},[1,2]);
figure,plot(vxy(:,:))

vper = ThreshCross(vxy(:,2),5,1);

figure();
plot(vxy(:,:));
Lines(vper(:,1),[],'g')
Lines(vper(:,2),[],'r')


pft = pfs_2d_theta(Trial,unitSubset);


unit = 103;
figure();
hold('on');
xlim([-500,500]);
ylim([-500,500]);
lgo = gobjects([0,1]);
rmap = plot(pft,unit,'mean','text');
gbins = cell([1,2]);
[gbins{:}] = meshgrid(pft.adata.bins{:});
plot(pft,unit,'mean','text');
contour(gbins{:},rmap',[1,1]);
for v = 1:size(vper,1),
    lgo(end+1) = plot(xyz(vper(v,1):vper(v,2),'head_back',1),   ...
                      xyz(vper(v,1):vper(v,2),'head_back',2),   ...
                      '.');
    waitforbuttonpress();
    if numel(lgo)>10,
        delete(lgo(1));
        lgo(1) = [];
    end
end


dist = [1:500]';
zoneBoundaries = [300,200,100,50];
zoneMat = bsxfun(@lt,dist,zoneBoundaries);
zoneId =  sum(zoneMat,2);
    
figure,plot(rhm.data);
Trial = MTATrial.validate('Ed01-20140709.cof.all');
xyz = Trial.load('xyz','trb');
rhm = fet_rhm(Trial);
lfp = Trial.load('lfp',2);
rhm.resample(lfp);

figure();
hold('on');
plot(nunity(rhm.data));
plot(nunity(lfp.data));


ang = create(MTADang,Trial,xyz);




 