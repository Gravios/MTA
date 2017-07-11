

Trial = MTATrial.validate('jg05-20120317.cof.all');

markers = {'spine_lower','pelvis_root','spine_middle','spine_upper',...
            'head_back',  'head_left',  'head_front',  'head_right'};



xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');

ang = create(MTADang,Trial,xyz);

vl = xyz.vel(markers,[1,2]);
vl.data(vl.data<1e-3) = 1e-4;
v = log10(abs(vl(:,1:8)));

edges = linspace(-3,2,64);
sbound = -30:30;
ixy = zeros([numel(sbound),size(v,2),size(v,2)]);

padding = [0,0];
padding = [-.5,.5];
vind = cast(resample(Trial.stc{'w'}+padding,xyz),'TimeSeries');
vind = logical(vind.data);
nind = sum(vind);


s = 1;
for m = 1:size(v,2)
    for o = 1:size(v,2)
        for shift = sbound
            [out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
            pxy = out./nind;
            px = histc(v(vind,m),xb);
            px = px(1:end-1)/nind;
            py = histc(circshift(v(vind,o),shift),yb);
            py = py(1:end-1)/nind;
            ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
            s = s+1;
        end
        s = 1;
    end
end




s = 1;
for m = 1:size(v,2)
    for o = 1:size(v,2)
        for shift = sbound
            [out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
            pxy = out./nind;
            px = histc(v(vind,m),xb);
            px = px(1:end-1)/nind;
            py = histc(circshift(v(vind,o),shift),yb);
            py = py(1:end-1)/nind;
            ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
            s = s+1;
        end
        s = 1;
    end
end



[mixy,sixy] = max(ixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);


figure,
subplot2(1,2,1,1);
imagesc(mixy(:,:));
caxis([0,2]);
colorbar
title('mutual information between marker speeds');
set(gca,'YtickMode','manual');
set(gca,'Ytick',1:8);
set(gca,'YtickLabelMode','manual');
set(gca,'YtickLabel',vl.model.ml('short'));
subplot2(1,2,1,2);
imagesc(sixy(:,:)/vl.sampleRate*1000);
colorbar
title('time lag of maximum mutual information (ms)')
set(gcf,'position',[520, 443, 1046, 289]);


stateLabels = {'all','rear','walk','turn','pause','groom','sit'};
states = 'arwnpms';

figure,hold on,
Trial.load('stc','hand_labeled_rev3_jg');

for s = states,
padding = [0,0];
padding = [-.5,.5];
vind = cast(resample(Trial.stc{s}+padding,xyz),'TimeSeries');
vind = logical(vind.data);
nind = sum(vind);

s = 1;
for m = 1,
for o = 7,
for shift = sbound
[out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
pxy = out./nind;
px = histc(v(vind,m),xb);
px = px(1:end-1)/nind;
py = histc(circshift(v(vind,o),shift),yb);
py = py(1:end-1)/nind;
ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
s = s+1;
end
s = 1;
end
end

plot(sbound./xyz.sampleRate,ixy(:,1,7))
end


legend(stateLabels)