function Trial = CorrectPointErrors(Trial,rb)

xyz = Trial.xyz(:,Trial.Model.gmi(rb.ml()),:);
errorPeriods = FindErrorPeriods(Trial);

edur =diff(errorPeriods,1,2);

%% Asymmetric
udur = unique(edur);
for eduration = 1:length(unique(edur)),
    pointErrors = errorPeriods(edur==udur(eduration),2);
    xyzSegs = GetSegs(xyz,pointErrors-udur(eduration),udur(eduration)+1);
    xyzSegs = reshape(xyzSegs,size(xyzSegs,1),size(xyzSegs,2),rb.N,3);
    marInd = zeros(udur(eduration),size(xyzSegs,2),rb.N);

    %if size(pointErrors,1)>10, continue;end


    for i = 2:udur(eduration)+1
        degMar = [];
        nonMar = [];
        tx =permute(reshape(repmat(xyzSegs(i-1:i,:,:,:),rb.N,1),2,rb.N,size(xyzSegs,2),rb.N,3),[1,3,2,4,5]);
        cx =permute(reshape(repmat(xyzSegs(i-1:i,:,:,:),rb.N,1),2,rb.N,size(xyzSegs,2),rb.N,3),[1,3,4,2,5]);
        cx(end,:,:,:,:) = tx(end,:,:,:,:);
        [~,marInd(i-1,:,:)] = min(sum(sum(diff(cx,1).^2,5),1),[],4);

        %% IF more than two markers are labeled with the same
        %% trajectory label as unfixable, marInd(...) = -ones;

        redMar = sq(sum(reshape(reshape(repmat([1:rb.N]',rb.N,size(xyzSegs,2)),rb.N,[])==repmat(reshape(sq(marInd(i-1,:,:))',[],1)',rb.N,1),rb.N,rb.N,[]),2));

        %        redMar = sq(sum(reshape(reshape(repmat([1:rb.N]',rb.N,size(mroo,2)),rb.N,[])==repmat(reshape(mroo,[],1)',rb.N,1),rb.N,rb.N,[]),2)) %% DELETE
        if length(find(redMar==1))<rb.N*size(xyzSegs,2),          %%DELETE
                                                                            %if length(find(redMar==1))<rb.N*size(mroo,2),             
            
            [nonMar nmi] = find(~redMar);
            [degMar dmi] = find(redMar==2);
            [~,imi] = find(redMar>2);
            %            degMarCmp = zeros(1,size(mroo,2));          %%DELETE
            degMarCmp = zeros(1,size(xyzSegs,2));
            degMarCmp(dmi) = degMar;
            degMarCmp(imi) = 0;

            %nonMarCmp = zeros(1,size(mroo,2));          %%DELETE
            nonMarCmp = zeros(1,size(xyzSegs,2));
            nonMarCmp(nmi) = nonMar;
            nonMarCmp(imi) = 0;

            %            degMarCmp(nonMarCmp==0)=0;

            rplMar = cumsum(repmat(degMarCmp,rb.N,1)==reshape(sq(marInd(i-1,:,:))',rb.N,[]));
            [rms,rmi] =sort(rplMar,1,'descend');
            rplMarInd = zeros(1,size(xyzSegs,2));
                                            %rplMarInd = zeros(1,size(mroo,2));                   %%DELETE
                                            %degMarInd = reshape((repmat(rplMarInd,rb.N,1)==repmat([1:rb.N]',1,size(mroo,2))).*repmat(nonMarCmp,rb.N,1),1,rb.N,size(mroo,2));
                                            %mmroo = reshape(mroo,1,rb.N,size(mroo,2));
                                            %mmroo(degMarInd~=0) = degMarInd(degMarInd~=0);
            rplMarInd(rms(1,:)==2) = rmi(1,rms(1,:)==2);
degMarInd = reshape(((repmat(rplMarInd,rb.N,1)==repmat([1:rb.N]',1,size(xyzSegs,2))).*repmat(nonMarCmp-degMarCmp,rb.N,1))',1,size(xyzSegs,2),rb.N);
%            degMarInd = reshape((repmat(rplMarInd,rb.N,1)==repmat([1:rb.N]',1,size(xyzSegs,2))).*repmat(nonMarCmp-degMarCmp,rb.N,1),1,size(xyzSegs,2),rb.N);
            %            degMarInd = reshape((repmat(rplMarInd,rb.N,1)==repmat([1:rb.N]',1,size(xyzSegs,2))).*repmat(nonMarCmp,rb.N,1),1,size(xyzSegs,2),rb.N);
            %            if udur(eduration)>1,keyboard;end
            %            keyboard
            marInd(i-1,:,:) = marInd(i-1,:,:)+degMarInd;
            if ~isempty(imi),marInd(i-1,imi,:) = reshape(repmat(1:rb.N,length(imi),1),1,length(imi),rb.N);end

        end

        
        %degMar = find(1<sum(repmat(sq(marInd(i-1,:,:)),rb.N,1)==repmat([1:rb.N]',size(xyzSegs,2),rb.N),3));
        %nonMar = find(1>sum(repmat(sq(marInd(i-1,:,:)),rb.N,1)==repmat([1:rb.N]',size(xyzSegs,2),rb.N),3));
        %% TODO - Prevent tajectory assignment to multiple markers
        for j = 1:size(xyzSegs,2)
            xyzSegs(i,j,[1:rb.N],:) = xyzSegs(i,j,marInd(i-1,j,:),:);
        end
    end
    xsegs = reshape(xyzSegs,[],rb.N,3);
    xind =[reshape(repmat(pointErrors'-udur(eduration),udur(eduration)+1,1),[],1)+repmat([0:udur(eduration)]',size(pointErrors,1),1)];
    Trial.xyz(xind,Trial.Model.gmi(rb.ml()),:) = xsegs;
end 

% $$$ 
% $$$ hbflr = Trial.transformOrigin('head_back','head_front',{'head_left','head_right'});
% $$$ hrlbf = Trial.transformOrigin('head_right','head_left',{'head_back','head_front'});
% $$$ sp3 = subplot(313);
% $$$ hold on
% $$$ Lines(errorPeriods(:,2),[-5,5],'r');
% $$$ plot([hbflr.transVec(:,:,2),hrlbf.transVec(:,:,2)])
% $$$ linkaxes([sp1,sp2,sp3],'xy')

