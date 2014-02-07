

pfr = MTAApfs(Trial,[],'rear',1,[],[30,30],[1.2,1.2],'xy');
pfw = MTAApfs(Trial,[],'walk',1,[],[30,30],[1.2,1.2],'xy');
pfh = MTAApfs(Trial,[],'hwalk',1,[],[30,30],[1.2,1.2],'xy');
pfl = MTAApfs(Trial,[],'lwalk',1,[],[30,30],[1.2,1.2],'xy');




% State height distribution
edges = [0:5:300];
Nrear = histc(Session.xyz(Session.stc{'r'},7,3),edges);

Nhwalk = histc(Session.xyz(Session.stc{'g'},7,3),edges);

Nlwalk = histc(Session.xyz(Session.stc{'l'},7,3),edges);

Ntheta = histc(Session.xyz(Session.stc{'t'},7,3),edges);


figure,
subplot2(4,1,1,1)
barh(edges,Nrear,'histc')
subplot2(4,1,2,1)
barh(edges,Nhwalk,'histc')
subplot2(4,1,3,1)
barh(edges,Nlwalk,'histc')
subplot2(4,1,4,1)
barh(edges,Ntheta,'histc')

