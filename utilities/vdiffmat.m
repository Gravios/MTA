
MTAstartup('cin','bach');
Trial = MTATrial('jg05-20120317');
Trial.load('xyz');
xy = Trial.xyz.data;

% xy: (frame(10),marker(4),xydim(2))
% $$$ xy = cat(3,...
% $$$   [209.5800  171.6553  122.2184   84.3474;...
% $$$    209.6121  171.7154  122.3711   84.4464;...
% $$$   209.6693  171.8336  122.4212   84.4685;...
% $$$   209.6453  171.8697  122.4939   84.5314;...
% $$$   209.7230  171.9547  122.4744   84.5605;...
% $$$   209.7480  171.9486  122.4922   84.5417;...
% $$$   209.7804  171.9544  122.4805   84.5451;...
% $$$   209.7508  171.9966  122.4549   84.5478;...
% $$$   209.7597  172.0115  122.4412   84.5347;...
% $$$   209.7683  172.0145  122.4501   84.5430],...
% $$$ ...
% $$$   [130.3834  118.3137  108.4700  109.2956;...
% $$$   130.3466  118.2474  108.2811  109.1514;...
% $$$   130.3495  118.1412  108.1969  109.1145;...
% $$$   130.3850  118.0852  108.1475  109.0753;...
% $$$   130.3660  117.9510  108.0802  109.0657;...
% $$$   130.3468  117.9083  108.0749  109.0886;...
% $$$   130.3301  117.8633  108.0513  109.1735;...
% $$$   130.3576  117.8145  108.0535  109.1944;...
% $$$   130.3088  117.7781  108.0529  109.1795;...
% $$$   130.2711  117.7605  108.0357  109.1685]);


xy = Trial.xyz.data;
% xy: (frame(10),marker(4),xydim(2))
nframe = size(xy,1); %nframe: number of frames (time)
nmar   =size(xy,2);  %nmar: number markers 
ndim   = size(xy,3); %ndim: number of spatial dimensions (xyz)
j =1:nmar;
tic
point_diffmat = sqrt(sum(cat(4,permute(reshape(repmat(xy(:,:,1)',nmar,1)-xy(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xy(:,:,2)',nmar,1)-xy(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xy(:,:,3)',nmar,1)-xy(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])).^2,4));
toc

ntf = 10.^[0:.5:5.5];

for bs = 1:1000,
for k = 1:numel(ntf),
    xy = Trial.xyz.data(1:ntf(k),:,:);
    j =1:nmar;
    nframe = size(xy,1); %nframe: number of frames (time)
    nmar   =size(xy,2);  %nmar: number markers 
    ndim   = size(xy,3); %ndim: number of spatial dimensions (xyz)


    %% Test For Loop with prealloc
    clear('point_diffmat')
    tic 
    point_diffmat = zeros([nframe,nmar,nmar,ndim]);
    for i=1:ndim,
        point_diffmat(:,:,:,i) = permute(reshape(repmat(xy(:,:,i)',nmar,1)-xy(:,j(ones(nmar,1),:),i).',[nmar,nmar,nframe]),[3,1,2]);
    end
    CompTime(k,1,bs) = toc;

    %% Test For Loop without prealloc
    clear('point_diffmat')
    tic 
    for i=1:ndim,
        point_diffmat(:,:,:,i) = permute(reshape(repmat(xy(:,:,i)',nmar,1)-xy(:,j(ones(nmar,1),:),i).',[nmar,nmar,nframe]),[3,1,2]);
    end
    CompTime(k,2,bs) = toc;


    %% Test double for loop for 
    clear('point_diffmat')
    tic
    point_diffmat = cat(4,permute(reshape(repmat(xy(:,:,1)',nmar,1)-xy(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xy(:,:,2)',nmar,1)-xy(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xy(:,:,3)',nmar,1)-xy(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2]));
    CompTime(k,3,bs) = toc;

    clear('point_diffmat')
    point_diffmat = zeros(xyz.size(1),xyz.size(2),xyz.size(2),3);
    for i=1:nmar,
        for j=1:nmar,
            point_diffmat(:,i,j,:) = xy(:,j,:)-xy(:,i,:);
        end
    end
    CompTime(k,4,bs) = toc;
    
end
end

figure,plot(ntf,CompTime(:,1),ntf,CompTime(:,2),ntf,CompTime(:,3))

point_diffmat = sqrt(sum(point_diffmat.^2,4));



nframe = size(xy,1); %nframe: number of frames (time)
nmar   =size(xy,2);  %nmar: number markers 
nmat   =size(xy,2)*size(xy,3);  %nmar: number markers 
ndim   = size(xy,3); %ndim: number of spatial dimensions (xyz)

txy = reshape(xy,[nframe,nmat]);
j =1:nmat;
point_diffmat = reshape(repmat(txy',nmat,1)-txy(:,j(ones(nmat,1),:)).',[nmat,nmat,nframe]);

point_diffmat = reshape(point_diffmat,[nframe,nmar,nmar,ndim*ndim]);
point_diffmat = sqrt(sum(point_diffmat.^2,5));
