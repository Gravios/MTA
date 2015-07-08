Trial = MTATrial('jg05-20120310');

units = select_units(Trial,25,'pyr');

xyz = Trial.load('xyz');

ufr = Trial.ufr.copy;

ufr.create(Trial,xyz,'theta',units,0.5);

[U,S,V] = svd(cov(ufr.data));

ucor = zeros(ufr.size);
for i= 1:ufr.size(2),
ucor(:,i) = sum(repmat(V(:,i)',[ufr.size(1),1]).*ufr.data,2);
end


tcor = ucor;
ucor = ufr.copy;
ucor.data = tcor;

uRcor = 3;
figure,imagesc(cov(ucor(Trial.stc{'r'},:)))

%[Ur,Sr,Vr] = svd(cov(ucor(Trial.stc{'r'},:)));

cov_state(:,:,1) = cov(ucor(Trial.stc{'r'},:));

for j = 1:ufr.size(2)
    mean_fet_state(:,j,1) = repmat(mean(ucor(nniz(ucor),j)),[ufr.size(1),1,1]);
end



mdl = fitlm(ucor,xyz(:,7,3));

ns = 1;
d_state = zeros(ucor.size(1),ns);
for i =  1:ns
    d_state(:,i) = -.5*log(det(cov_state(:,:,i)))...
        -.5*dot(((ucor.data(:,:)...
        -mean_fet_state(:,:,i))/cov_state(:,:,i))',(ucor.data(:,:)...
        -mean_fet_state(:,:,i))');
end

