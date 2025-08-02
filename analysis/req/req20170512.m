
Trial = MTATrial.validate('jg05-20120317.cof.all');
stc =  Trial.load('stc','hand_labeled_rev3_jg');
fet = fet_mis(Trial,Trial.xyz.sampleRate);


ufet = fet.copy;
for i = 1:size(fet,2),
    [ufet.data(nniz(fet),i)] = make_uniform_unit_distribution(fet(nniz(fet),i));
end

fufet = ufet.copy;
fufet.filter('ButFilter',3,2.4,'low');

states = {'walk','turn','pause','rear'};
sclr = 'bgcr';
sp = [];
figure,
sp(end+1)=subplot2(10,1,1:5,1);
imagesc([1:size(fet,1)]./fet.sampleRate,1:18,ufet.data');
sp(end+1)=subplot2(10,1,[6:7],1);
plot([1:size(fet,1)-1]./fet.sampleRate,sqrt(sum(diff(ufet.data).^2,2)))
sp(end+1)=subplot2(10,1,[8:10],1);
plotSTC(stc,1,'text',states,sclr,[],false);

linkaxes(sp,'x');


cov(ufet(stc{'r'},:))
mean(ufet(stc{'r'},:))
%mean(ufet(stc{'w'},:))
%mean(ufet(stc{'n'},:))


mean(ufet(stc{'n'},:))

stsMeanDistMat = [];
sts = 'wrnpms';
nsts = numel(sts);
for o = 1:nsts,
    for t = 1:nsts,
        stsMeanDistMat(o,t) = sqrt(sum([mean(ufet(stc{sts(o)},:))-mean(ufet(stc{sts(t)},:))].^2));
    end
end

        
        
        
        
        
        