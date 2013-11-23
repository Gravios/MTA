function Trial = CorrectRigidBody(Trial,rb,pose)

%[rb,imd_fit,~,model_index] = RBErrorModel(Trial,markerSet);
xyz = Trial.xyz(:,Trial.Model.gmi(rb.ml()),:);
dist = Trial.ang(:,Trial.Model.gmi(rb.ml()),Trial.Model.gmi(rb.ml()),3);
depth = size(rb.imdMean,3);

%% Get all Possible Marker Swaps
error_key = perms([1:rb.N]);
[~,error_remap] = sort(error_key,2);


%% (Model Fit Score) Used to find the frames which need repair
imd_fit = zeros(size(Trial.xyz,1),depth);
dist_fit = zeros(size(Trial.xyz,1),depth);
for i = 1:size(xyz,1),
    dist_fit(i,:) = sq(sum(sum(reshape(repmat(abs(sq(dist(i,:,:))-rb.imdMean(:,:,depth)),1,depth),rb.N,rb.N,depth)-rb.imdStd(:,:,:))));
end
imd_fit = sort(abs(dist_fit),2,'ascend');
%%%


model_swap_perms = reshape(pose(error_key,:),size(error_key,1),rb.N,3);

mimo = imo(model_swap_perms,rb);
aimo = imo(xyz,rb);

error_chk = zeros(size(mimo,1),size(error_key,1));
error_fet = zeros(size(xyz,1),size(error_key,1));
error_map = zeros(size(error_key,1),1);
for i = 1:size(error_key,1),
    error_chk(:,i) = mimo(:,error_key(i,1),error_key(i,2),error_key(i,3),error_key(i,4),1);
    error_fet(:,i) = aimo(:,error_key(i,1),error_key(i,2),error_key(i,3),error_key(i,4),1);
    error_map(i) = find(ismember(error_key,error_key(i,:),'rows'));
end

[~,error_chk_ind] = sort(error_chk,2);
[~,error_fet_ind] = sort(error_fet,2);
eci = permute(reshape(repmat(error_chk_ind,1,size(xyz,1)),size(error_chk_ind,1),size(error_chk_ind,2),[]),[3,1,2]);
efi = permute(reshape(repmat(error_fet_ind,1,size(error_chk_ind,1)),size(xyz,1),size(error_fet_ind,2),[]),[1 3 2]);

error_fit = sum(eci==efi,3);
[~,error_typ] = sort(error_fit,2,'descend');
error_ind = find(error_typ(:,1)~=find(ismember(error_key,[1:rb.N],'rows')));
for i = 1:size(error_ind,1),
    xyz(error_ind(i),:,:) = xyz(error_ind(i),error_remap(error_typ(error_ind(i),1),:),:);
end


TempTrial = Trial;
TempTrial.xyz(:,Trial.Model.gmi(rb.ml()),:) = xyz;
TempTrial = TempTrial.load_ang(1);
TempTrial.updateModel();

%% (Model Fit Score) Used to find the frames which need repair
trb = TempTrial.Model.rb(rb.ml());
temp_imd_fit = zeros(size(xyz,1),depth);
dist = Trial.ang(:,TempTrial.Model.gmi(rb.ml()),TempTrial.Model.gmi(rb.ml()),3);
for i = 1:size(xyz,1),
    dist_fit(i,:) = sq(sum(sum(reshape(repmat(abs(sq(dist(i,:,:))-rb.imdMean(:,:,depth)),1,depth),rb.N,rb.N,depth)-rb.imdStd(:,:,:))));
end
temp_imd_fit = sort(abs(dist_fit),2,'ascend');


improved_imd_fit = zeros(size(xyz,1),1);
improved_imd_fit(temp_imd_fit(:,depth)<=imd_fit(:,depth)) = 1;
improved_imd_fit(find(error_typ(:,1)==find(ismember(error_key,[1:rb.N],'rows'))))=0;

Trial.xyz(improved_imd_fit==1,Trial.Model.gmi(rb.ml()),:) = xyz(improved_imd_fit==1,:,:);
Trial = Trial.load_ang(1);





