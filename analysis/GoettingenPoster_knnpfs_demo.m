
Trial = MTATrial('jg05-20120310');


if Trial.xyz.isempty, Trial.xyz.load(Trial);end
sstpos = sq(Trial.xyz(Trial.stc{'w'},7,[1,2])); 
Trial.ufr.create(Trial,Trial.xyz,'walk',10,0.2);
sstufr = Trial.ufr(Trial.stc{'w'},:);

rates = unique(sstufr);
nrates = numel(rates);

for i=rates'
    sstufr(sstufr==i) = find(rates==i);
end
clrs = jet(2*nrates);
clrs = clrs([1,nrates-1:end],:);

figure
scatter(sstpos(:,1),sstpos(:,2),10,clrs(sstufr,:),'filled');