

Trial = MTATrial('jg05-20120317');

units = [];

pfx = MTAApfs(Trial,units,'theta',1,[],[30,30],[1.2,1.2],'xy');
pfs = MTAApfs(Trial,units,'theta',1,[],[30,30,30],[1.2,1.2,1.2],'xyz');
pfr = MTAApfs(Trial,units,'theta',1,[],[30,30],[1.2,1.2],'pfcrz');

units = pfs.data.clu;

Trial.spk.create(Trial,Trial.xyz.sampleRate,'theta');
%Trial.ufr.create(Trial,pfs,'theta');
Trial.ufr = MTADufr([]);
Trial.ufr.create(Trial,Trial.xyz,'theta',units);




u =29;

uind = u==units;
ufrmax = max(Trial.ufr(:,uind));
uthr = [.6,1].*ufrmax;
ufrind = Trial.ufr(:,uind)>uthr(1)&Trial.ufr(:,uind)<uthr(2);
c = jet(5);
ufrvals = Trial.ufr(ufrind,uind);
ufrc = ceil((ufrvals-uthr(1)).*5/abs(diff(uthr)));

ufrindex = find(ufrind);
iind = randi(numel(ufrindex),[10000,1]);
jind = randi(numel(ufrindex),[10000,1]);
dmat = sqrt(sum((Trial.xyz(ufrindex(iind),7,:)-Trial.xyz(ufrindex(jind),7,:)).^2,3));

%figure,
subplot(151),cla
pfs.plot(u)
zlim([0,500])
subplot(152),cla
pfx.plot(u)
subplot(153),cla
pfr.plot(u)
subplot(154),cla
% figure
scatter3(Trial.xyz(ufrind,7,1),Trial.xyz(ufrind,7,2),Trial.xyz(ufrind,7,3),20,c(ufrc,:));%,'MarkerFace','filled')
xlim([-500,500]),ylim([-500,500]),zlim([0,500])
subplot(155),cla
hist(dmat,1000),Lines(mean(dmat),[],'k');
kurtosis(dmat)
skewness(dmat)


figure
xk = [];
xs = [];
    u = 95
    for i=1:10,
    uind = u==units;
    ufrmax = max(Trial.ufr(:,uind));
    uthr = [.08*i,1].*ufrmax;
    ufrind =2 Trial.ufr(:,uind)>uthr(1)&Trial.ufr(:,uind)<uthr(2);
    c = jet(5);
    ufrvals = Trial.ufr(ufrind,uind);
    ufrc = ceil((ufrvals-uthr(1)).*5/abs(diff(uthr)));
    %subplot(10,1,i),
    % edges = [0,25,50,75,100,125,150,175,200,225,250,275,300];
    % N = histc(Trial.xyz(ufrind,7,3),edges);
    % bar(edges,N,'histc')
    %hist2([ufrvals,Trial.xyz(ufrind,7,3)],50,50);
    xk(i,:) = kurtosis(sq(Trial.xyz(ufrind,7,:)));
    xs(i,:) = skewness(sq(Trial.xyz(ufrind,7,:)));
    end
%plot3(Trial.xyz(ufrind,7,1),Trial.xyz(ufrind,7,2),Trial.xyz(ufrind,7,3),'.')
%plot3(Trial.xyz(Trial.spk(u),7,1),Trial.xyz(Trial.spk(u),7,2),Trial.xyz(Trial.spk(u),7,3),'.')


