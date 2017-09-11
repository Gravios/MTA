function Trial = CorrectPointErrors(Trial,rb)

if Trial.xyz.isempty,Trial.xyz.load(Trial);end

xyz = Trial.xyz(:,rb.ml,:);

errorPeriods = FindErrorPeriods(Trial);

edur =diff(errorPeriods,1,2);

%% Asymmetric
udur = unique(edur);
for eduration = 1:length(unique(edur)),
    pointErrors = errorPeriods(edur==udur(eduration),2);
    xyzSegs = GetSegs(xyz,pointErrors-udur(eduration),udur(eduration)+1);
    xyzSegs = reshape(xyzSegs,size(xyzSegs,1),size(xyzSegs,2),rb.N,3);
    marInd = zeros(udur(eduration),size(xyzSegs,2),rb.N);

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

        if length(find(redMar==1))<rb.N*size(xyzSegs,2),          

            
            [nonMar nmi] = find(~redMar);
            [degMar dmi] = find(redMar==2);
            [~,imi] = find(redMar>2);

            degMarCmp = zeros(1,size(xyzSegs,2));
            degMarCmp(dmi) = degMar;
            degMarCmp(imi) = 0;

            nonMarCmp = zeros(1,size(xyzSegs,2));
            nonMarCmp(nmi) = nonMar;
            nonMarCmp(imi) = 0;


            rplMar = cumsum(repmat(degMarCmp,rb.N,1)==reshape(sq(marInd(i-1,:,:))',rb.N,[]));
            [rms,rmi] =sort(rplMar,1,'descend');
            rplMarInd = zeros(1,size(xyzSegs,2));

            rplMarInd(rms(1,:)==2) = rmi(1,rms(1,:)==2);
            degMarInd = reshape(((repmat(rplMarInd,rb.N,1)==repmat([1:rb.N]',1,size(xyzSegs,2))).*repmat(nonMarCmp-degMarCmp,rb.N,1))',1,size(xyzSegs,2),rb.N);

            marInd(i-1,:,:) = marInd(i-1,:,:)+degMarInd;
            if ~isempty(imi),marInd(i-1,imi,:) = reshape(repmat(1:rb.N,length(imi),1),1,length(imi),rb.N);end

        end

        %% TODO - Prevent tajectory assignment to multiple markers
        for j = 1:size(xyzSegs,2)
            xyzSegs(i,j,[1:rb.N],:) = xyzSegs(i,j,marInd(i-1,j,:),:);
        end
    end
    xsegs = reshape(xyzSegs,[],rb.N,3);
    xind =[reshape(repmat(pointErrors'-udur(eduration),udur(eduration)+1,1),[],1)+repmat([0:udur(eduration)]',size(pointErrors,1),1)];
    Trial.xyz(xind,rb.ml,:) = xsegs;
end 


