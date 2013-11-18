
t = MTATrial('jg05-20120317');
t.ang.load(t.sync);


windows = [9,27,63,81,121];
dims = [1,2];
dsf = 4;

sl = t.xyz.size(1);
sl = sl-mod(sl,dsf);

cons = [1,2,3,4,5;...
        2,3,4,5,7];
    
fdda = zeros(sl/dsf,size(cons,2),numel(windows),numel(dims));
for n = 1:numel(windows)
    tfdda = zeros(windows(n),sl/dsf,size(cons,2),numel(dims));
    for c =cons
        tfdda(:,:,c(1)==cons(1,:),:) = GetSegs(sq(t.ang(1:sl,c(1),c(2),dims)),1:dsf:sl,windows(n),0);
    end
    fdda(:,:,n,:) = sq(median(cumsum(circ_dist(tfdda,cat(1,tfdda(1,:,:,:),tfdda(1:end-1,:,:,:))))));
end
mfdda = MTADxyz([],[],fdda,t.ang.sampleRate/dsf);





mxyz = MTADxyz([],[],t.xyz(round(linspace(1,t.xyz.size(1),mfdda.size(1))),:,:),t.xyz.sampleRate/dsf);

mang = MTADxyz([],[],t.ang(round(linspace(1,t.ang.size(1),mfdda.size(1))),:,:,:),t.xyz.sampleRate/dsf);

vel = t.vel;
mvel = MTADxyz([],[],vel(round(linspace(1,t.xyz.size(1)-1,mfdda.size(1))),:),t.xyz.sampleRate/dsf);

imagesc(log10(mvel(:,:)+0.01)')

mv = mean(mvel(mvel(:,1)~=0,:));
sv = std(mvel(mvel(:,1)~=0,:));
uv = (mvel(:,:)-repmat(mv,mvel.size(1),1))./repmat(sv,mvel.size(1),1);

imagesc(uv'),caxis([0,3])

mz = mean(mxyz(mxyz(:,1,1)~=0,:,3));
sz = std(mxyz(mxyz(:,1,1)~=0,:,3));
uz = (mxyz(:,:,3)-repmat(mz,mxyz.size(1),1))./repmat(sz,mxyz.size(1),1);

imagesc(uz'),caxis([0,3])


sp=[];
figure,
sp(1)=subplot2(10,1,1:9,1);imagesc(linspace(1,t.ang.size(1),mfdda.size(1))./t.ang.sampleRate,1:28,cat(2,2.^uv-1,uz,2.^abs(sq(mfdda(:,:,3,2)))-1,2.^abs(sq(mfdda(:,:,3,1)))-1)'),caxis([-.5,3])
sp(2)=subplot2(10,1,10,1); image(linspace(1,t.ang.size(1),mfdda.size(1))./t.ang.sampleRate,1,stb)
linkaxes(sp,'x')


stb = ones(1,mxyz.size(1),3);
color = [0,1,0;...
         1,0,0;...
         0,0,1];
states = {'b','r','w'};      
for i = 1:3
sts = t.stc{states{i},mxyz.sampleRate}.data';
for s = sts,
    stb(1,s(1):s(2),:) = repmat(color(i,:),diff(s)+1,1);
end
end
    


