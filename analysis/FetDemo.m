
s = MTASession('jg05-20120317');
t = MTATrial(s);



plot(circ_dist(t.ang(t.stc{'w'}(:,:),1,2,1),t.ang(t.stc{'w'}(:,:),3,4,1))),ylim([-pi,pi])
hold on,
Lines(cumsum([0;diff(t.stc{'w'}(:,:),1,2)+1]),[],'k')
plot(circ_dist(t.ang(t.stc{'w'}(:,:),1,2,1),t.ang(t.stc{'w'}(:,:),5,7,1)),'r'),ylim([-pi,pi])
plot(circ_dist(t.ang(t.stc{'w'}(:,:),3,4,1),t.ang(t.stc{'w'}(:,:),5,7,1)),'c'),ylim([-pi,pi])
plot(circ_dist(t.ang(t.stc{'w'}(:,:),3,4,1),t.ang(t.stc{'w'}(:,:),4,5,1)),'m'),ylim([-pi,pi])


dper = diff(t.stc{'w'}(:,:),1,2);
ba = zeros(sum(dper+1),5);
for i = 1:5,
ba(:,i) = circ_dist(t.ang(t.stc{'w'}(:,:),1,2,1),t.ang(t.stc{'w'}(:,:),1,i+2,1));
end
figure,imagesc(ba')

circ_dist(t.ang(t.stc{'w'}(:,:),1,2,1),t.ang(t.stc{'w'}(:,:),3,4,1))),ylim([-pi,pi])
hold on,
Lines(cumsum([0;dper+1]),[],'k')
plot(circ_dist(t.ang(t.stc{'w'}(:,:),1,2,1),t.ang(t.stc{'w'}(:,:),5,7,1)),'r'),ylim([-pi,pi])
plot(circ_dist(t.ang(t.stc{'w'}(:,:),3,4,1),t.ang(t.stc{'w'}(:,:),5,7,1)),'c'),ylim([-pi,pi])
plot(circ_dist(t.ang(t.stc{'w'}(:,:),3,4,1),t.ang(t.stc{'w'}(:,:),4,5,1)),'m'),ylim([-pi,pi])



figure,
plot(sq(log10(mdx(:,[1,3,7],3)))),
Lines(t.stc{'w',mdx.sampleRate}(:),[],'k');
Lines(t.stc{'r',mdx.sampleRate}(:),[],'r');

