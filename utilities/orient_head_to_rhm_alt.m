function orient_head_to_rhm_alt(Trial,xyz,ncpChannel)


Trial = 'Ed01-20140709.cof.all';
Trial = 'Ed01-20140707.cof.all';

Trial = 'Ed03-20140624.cof.all';
Trial = 'Ed03-20140625.cof.all';

Trial = 'Ed05-20140528.cof.all';
Trial = 'Ed05-20140529.ont.all';

Trial = 'Ed10-20140816.cof.all';
Trial = 'Ed10-20140817.cof.gnd';

Trials = {'Ed03-20140624.cof.all','Ed03-20140625.cof.all'};

for t = 1:numel(Trials);
Trial = MTATrial.validate(Trials{t});
ncpChannel = [];
%ncpChannel = 66;
xyz = preproc_xyz(Trial,'trb');
ncp = fet_ncp(Trial,[],[],ncpChannel);

vxy = xyz.vel({'hcom','head_front'},[1,2]);
vxy.data(vxy(:)<1e-3) = 1e-3;
vxy = log10(vxy.data);


% Need basis vectors from transtform_rigidbody.m
blen = 45;

xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[3]./(xyz.sampleRate/2),'low'));
xyz.data(~nniz(xyz(:,1,1)),end,:) = 0;


% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
nm = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',  [0.5,1,0.5],[],nm);
nm = nz.*45+xyz(:,'hcom',:);
xyz.addMarker('hrhm',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(xyz(:,'htx',:)-xyz(:,'hcom',:),xyz(:,'head_back',:)-xyz(:,'hcom',:));
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+xyz(:,'hcom',:);
xyz.addMarker('hrx',  [0.5,1,0.5],[],nm);


% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(xyz(:,'hrx',:)-xyz(:,'hcom',:),xyz(:,'htx',:)-xyz(:,'hcom',:));
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+xyz(:,'hcom',:);
xyz.addMarker('hbx',  [0.5,1,0.5],[],nm);

evec = bsxfun(@minus,xyz(:,{'hbx','hrx','htx'},:),xyz(:,{'hcom'},:));


pind = nniz(evec)&vxy(:,2)>0;

dispFigs = true;
moffset = 1000;

markers = {'hcom','head_front'};

% Find Best orientation (good luck)

%% Display the Original basis for the head
if dispFigs,
    hfig = figure(20392039);
    set(hfig,'units','centimeters')
    set(hfig,'Position',[11 0 32 25])
    set(hfig,'PaperPositionMode','auto');
    clf;
    subplot(2,2,1);
    hold('on');
    sind = find(pind,1,'first')+moffset;    
    hline = gco;
    markers = {'head_back','head_left','head_front','head_right','hbx','hrx','htx'};
    mcolors = 'bgmrbrcb';
    mlineWidth = ones([1,numel(mcolors)]); mlineWidth(end)=3;
    mlineMarker = repmat({'none'},[1,numel(mcolors)]);
    mlineMarker{end} = 'o';
    for m = 1:numel(markers)
        hline(end+1) = line(xyz(sind,{'hcom',markers{m}},1),... X
                            xyz(sind,{'hcom',markers{m}},2),... Y
                            xyz(sind,{'hcom',markers{m}},3),... Z
                            'Color',mcolors(m),...
                            'LineWidth',mlineWidth(m),...
                            'Marker',mlineMarker{m});
    end
    
    title('Vector Search Space');
    daspect([1,1,1]);
end

ind = vxy(:,2)>0.3;
aind = Trial.stc{'a'};
aind.cast('TimeSeries');
aind.resample(xyz);
aind.data(isnan(aind.data))=0;
ind = logical(aind.data)&ind;

ind = MTADepoch([],[],ind,xyz.sampleRate,xyz.sync,xyz.origin,'TimeSeries');


%% Rotate a basis around the unit sphere
spw = [];
cpw = [];
rpw = [];
npw = [];


rots = 2.5:5:180;
tots = 2.5:5:360;

% $$$ rots = 15:15:180;
% $$$ tots = 0:15:360;
% $$$ 
% $$$ rots = 5:40:180;
% $$$ tots = 5:40:360;
% $$$ 

for d =  1:numel(rots)
% SET rotation axis to z-axis (htx)
% ROTATE around z-axis (htx)
    m = 1;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(evec(:,m,ax_ord)),sqrt(sum(evec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
% CREATE rotation matrix
    j = 1:3;
    rot_ang = deg2rad(rots(d));
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

% SET matrix
    ovec = evec(:,2,:);
% CREATE new basis vector from rotation
    nvec = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,2,1]);
    nvec = cat(2,evec,nvec);

% SET rotation axis to new vector
    m = 2;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(nvec(:,m,ax_ord)),sqrt(sum(nvec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
         -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);

    for t =  1:numel(tots)
% ROTATE sep wise around new vector
        j =1:3;  
        rot_ang = deg2rad(tots(t));
        head_rotMat =  cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
                      +sin(rot_ang)*head_cpm...
                      +(1-cos(rot_ang))*head_kron;

        % Rotated marker;???
        ovec = nvec(:,4,:);
        ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;

        % Create new marker
        tmark = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]), ...
                                                 [2,3,1]),2),[3,1,2])+sq(xyz(:,'hcom',:));

        sxyz = MTADfet.encapsulate(Trial,cat(2,xyz(:,'fhcom',:),permute(tmark,[1,3,2])),xyz.sampleRate,...
                                   'sxyz','sxyz','s');
        sxyz.filter('RectFilter',3,4);
        sxyz.data = diff([0;diff(sqrt(sum(diff(sxyz.data,[],2).^2,3)));0]);        
        rxyz = sxyz.copy();
        rxyz.data(~nniz(xyz(:,1,1))) = 0;
        rxyz.filter('ButFilter',3,[6,14],'bandpass');
        xsegs = nansum(GetSegs(rxyz(pind),[1:round(xyz.sampleRate/2):sum(pind)]-round(xyz.sampleRate),...
                                         round(xyz.sampleRate)).^2);
        spw(d,t) = nanmedian(xsegs);

        [ys,fs,ts] = mtcsdglong([sxyz.data,ncp.data],2^8,sxyz.sampleRate,2^7,2^6,[],[],[],[5,14]);
        if d ==1&&t==1,
            ts = ts+(2^6/2)/sxyz.sampleRate;
            ssr = 1/diff(ts(1:2));
            pad = round([ts(1),size(sxyz,1)./sxyz.sampleRate-ts(end)].*ssr)-[1,0];
            szy = size(ys);

            resample(ind,ssr);
            if size(ys,1)>size(ind,1),
                ind.data = cat(1,ind.data,zeros([size(ys,1)-size(ind,1),1]));
            elseif size(ys,1)<size(ind,1),
                ind.data = ind.data(1:size(ys,1));
            end
            ind = ind.data;

        end 
        yss = sq(mean(cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),2));
        cpw(d,t) = nanmean(abs(yss(ind,1,2)))./nanmean(sqrt(yss(ind,1,1).*yss(ind,2,2)));
        rpw(d,t) = log10(nanmean(abs(yss(ind,1,1))));
        if dispFigs,
            plot3(tmark(sind,1),tmark(sind,2),tmark(sind,3),'*r');
            drawnow();
        end
    end%for t
end%for d

save(['/storage/gravio/data/project/general/analysis/orient_head_to_rhm_alt-',Trial.filebase,'.mat'],...
     'rots','tots','spw','cpw','rpw');
end
Trial = 'Ed05-20140528.cof.all';
Trial = 'Ed05-20140529.ont.all';

Trial = 'Ed01-20140709.cof.all';
Trial = 'Ed01-20140707.cof.all';

Trial = 'Ed03-20140624.cof.all';
Trial = 'Ed03-20140625.cof.all';


Trial = 'Ed10-20140816.cof.all';
Trial = 'Ed10-20140817.cof.gnd';
Trial = MTATrial.validate(Trial);
load(['/storage/gravio/data/project/general/analysis/orient_head_to_rhm_alt-',Trial.filebase,'.mat']);


figure,
subplot(221);imagesc(rots,tots,spw');%caxis([0.07,0.15]); estimate
subplot(222);imagesc(rots,tots,cpw'),%caxis([0.4,0.55]);
subplot(223);imagesc(rots,tots,rpw'),%caxis([0.4,0.55]);
subplot(224);plot([tots;tots]',[spw(19,:);cpw(19,:)]')






%subplot(224);imagesc(rots,tots,npw'),%caxis([0.4,0.55]);
% SMOOTH search space

% $$$ sspw = [flipud(circshift(spw',round(size(spw,2)/2))),fliplr(spw')];
% $$$ sspw = [sspw,circshift(sspw,round(size(sspw,2)/8))];



%% Automatic orientation detection
% spw is the map where the peaks represents the front and back of
% the head. The indicies (d,t) are the two angle thingies you need 
% to get the final orientation vector. 

Trials = {'Ed05-20140528.cof.all','Ed05-20140529.ont.all',...
          'Ed01-20140709.cof.all','Ed01-20140707.cof.all',...
          'Ed03-20140624.cof.all','Ed03-20140625.cof.all',...
          'Ed10-20140816.cof.all','Ed10-20140817.cof.gnd'};

for t = 1:numel(Trials);
    Trial = MTATrial.validate(Trials{t});
    load(['/storage/gravio/data/project/general/analysis/orient_head_to_rhm_alt-',Trial.filebase,'.mat']);

srpw = repmat(spw,[3,3]);
sspw = repmat(cpw,[3,3]);

%figure,imagesc(sspw')

SmoothingWeights = [3,3];
msize = size(sspw);
ndims = numel(msize);

sind = cell(1,2);
for i = 1:ndims,
    sind{i} = linspace(-round(msize(i)/2),round(msize(i)/2),msize(i));
end
[sind{:}] = ndgrid(sind{:});
for i = 1:ndims,
    sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));

srpw = convn(srpw,Smoother,'same');
sspw = convn(sspw,Smoother,'same');

frpw = srpw(size(spw,1)+1:size(spw,1)*2,size(spw,2)+1:size(spw,2)*2);
fspw = sspw(size(spw,1)+1:size(spw,1)*2,size(spw,2)+1:size(spw,2)*2);

% $$$ if dispFigs,
% $$$     subplot(222);
% $$$     imagesc(rots,tots,fspw');
% $$$     title('Rhythmic(6-13Hz) mean PSD')
% $$$     xlabel('\phi (rad)')
% $$$     ylabel('\theta (rad)');
% $$$ end
[mins,minv] = LocalMinimaN(-sspw,-nanmean(sspw(:)),10);
[minsr,minvr] = LocalMinimaN(-rspw,-nanmean(rspw(:)),10);

mins = bsxfun(@minus,mins,[numel(rots),numel(tots)]);


minv(~[mins(:,1)>0 & mins(:,2)>0 & mins(:,1)<numel(rots) & mins(:,2)<numel(tots)]) = [];
mins(~[mins(:,1)>0 & mins(:,2)>0 & mins(:,1)<numel(rots) & mins(:,2)<numel(tots)],:) = [];

minvr(~[minsr(:,1)>0 & minsr(:,2)>0 & minsr(:,1)<numel(rots) & minsr(:,2)<numel(tots)]) = [];
minsr(~[minsr(:,1)>0 & minsr(:,2)>0 & minsr(:,1)<numel(rots) & minsr(:,2)<numel(tots)],:) = [];

% Get best vector set
cang(t) = tots(mins(end));
rncd(t) = circ_dist(tots(mins(end)),tots(minsr(end)));

end

bvecs = [];

for dt = mins'
% ROTATE around z-axis (htx)
    m = 1;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(evec(:,m,ax_ord)),sqrt(sum(evec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);

    j = 1:3;
    rot_ang = deg2rad(rots(dt(1)));    
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

    % Rotated marker; Head back
    ovec = evec(:,2,:);

    nvec = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,2,1]);
    nvec = cat(2,evec,nvec);

    m = 2;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(nvec(:,m,ax_ord)),sqrt(sum(nvec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
    
    j =1:3;  
    rot_ang = deg2rad(tots(dt(2)));
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

    % Rotated marker;???
    ovec = nvec(:,4,:);    
    ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;

    % Create new marker
    tmark = permute(permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),...
                                                      [size(head_norm,1),3,3]), ...
                                             [2,3,1]),2),...
                            [3,1,2])+sq(xyz(:,'hcom',:)),...
                    [1,3,2]);
    bvecs = cat(2,bvecs,tmark);

end



% $$$ bvecs = [];
% $$$ 
% $$$ for dt = mins'
% $$$     m = 3;
% $$$     ax_ord = [1,2,3];
% $$$     j =1:3;
% $$$     head_norm = bsxfun(@rdivide,sq(evec(ind,m,ax_ord)),sqrt(sum(evec(ind,m,ax_ord).^2,3)));
% $$$     head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
% $$$     j = [ 0,-1, 1;...
% $$$           1, 0,-1;...
% $$$          -1, 1, 0];
% $$$     k = [1,3,2;...
% $$$          3,1,1;...
% $$$          2,1,1];
% $$$     head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
% $$$ 
% $$$     j = 1:3;  
% $$$     rot_ang = deg2rad(rots(dt(1)));
% $$$     head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
% $$$                   +sin(rot_ang)*head_cpm...
% $$$                   +(1-cos(rot_ang))*head_kron;
% $$$ 
% $$$     % Rotated marker; Head Back
% $$$     ovec = evec(ind,1,:);
% $$$ 
% $$$     nvec = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,2,1]);
% $$$     nvec = cat(2,evec(ind,:,:),nvec);
% $$$ 
% $$$     m = 4;
% $$$     ax_ord = [1,2,3];
% $$$     j =1:3;
% $$$     head_norm = bsxfun(@rdivide,sq(nvec(:,m,ax_ord)),sqrt(sum(nvec(:,m,ax_ord).^2,3)));
% $$$     head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
% $$$     j = [ 0,-1, 1;...
% $$$           1, 0,-1;...
% $$$          -1, 1, 0];
% $$$     k = [1,3,2;...
% $$$          3,1,1;...
% $$$          2,1,1];
% $$$     head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
% $$$ 
% $$$     j =1:3;  
% $$$     rot_ang = deg2rad(tots(dt(2)));
% $$$     head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
% $$$         +sin(rot_ang)*head_cpm...
% $$$         +(1-cos(rot_ang))*head_kron;
% $$$ 
% $$$     % Rotated marker;
% $$$     %ovec = nvec(ind,2,:);         head right ???
% $$$     ovec = nvec(:,2,:);
% $$$     ovec = cross(nvec(:,2,:),nvec(:,4,:));
% $$$     ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;
% $$$ 
% $$$     tmark = permute(permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2])+sq(xyz(ind,'hcom',:)),[1,3,2]);
% $$$     
% $$$     bvecs = cat(2,bvecs,tmark);
% $$$ end




dout = bvecs(:,min(abs(minv))==abs(minv),:);
ovecs = bvecs;
% $$$ if dispFigs,
% $$$     suptitle([filename,' SessionId: ',num2str(ses)])
% $$$     print(hfig,'-depsc2',fullfile(odir,ofilename,[ofilename,'-sessionId-',num2str(ses),'.eps']));
% $$$     print(hfig,'-dpng',fullfile(odir,ofilename,[ofilename,'-sessionId-',num2str(ses),'.png']));
% $$$     delete(hfig)
% $$$ end


%% Save the file


xyz.addMarker('hfront',  [0.5,1,0.5],[],zeros(size(xyz,1),1,3));
xyz.data(:,end,:) = dout;

nxyz = Trial.load('xyz','trb');
nxyz.addMarker('hfront',  [0.5,1,0.5],[],zeros(size(xyz,1),1,3));
nxyz.data(ind,end,:) = dout;

ang = create(MTADang,Trial,xyz);

figure,plot(circ_dist(ang(:,'head_back','head_front',1),ang(:,'head_back','hfront',1)))
figure,plot(circ_dist(ang(:,'head_back','head_front',2),ang(:,'head_back','hfront',2)))
figure,plot(ang(:,'head_front','hfront',3))

%% Display the Original basis for the head
if dispFigs,
    hfig = figure(20392010);
    set(hfig,'units','centimeters')
    set(hfig,'Position',[11 0 32 25])
    set(hfig,'PaperPositionMode','auto');
    clf();
% $$$     subplot(2,2,1);
    hold('on');
    sind = 98000;
    hline = gco;
    markers = {'head_back','head_left','head_front','head_right','hbx','hrx','htx','hfront'},
    mcolors = 'bgmrbrcb';
    mlineWidth = ones([1,numel(mcolors)]); mlineWidth(end)=3;
    mlineMarker = repmat({'none'},[1,numel(mcolors)]);
    mlineMarker{end} = 'o';
    for m = 1:numel(markers)
        hline(end+1) = line(xyz(sind,{'hcom',markers{m}},1),... X
                            xyz(sind,{'hcom',markers{m}},2),... Y
                            xyz(sind,{'hcom',markers{m}},3),... Z
                            'Color',mcolors(m),...
                            'LineWidth',mlineWidth(m),...
                            'Marker',mlineMarker{m});
    end
    for m = 1:size(bvecs,2)
        line([xyz(sind,{'hcom'},1),ovecs(sind,m,1)],... X
                            [xyz(sind,{'hcom'},2),ovecs(sind,m,2)],... Y
                            [xyz(sind,{'hcom'},3),ovecs(sind,m,3)],... Z
                            'Color','k',...mcolors(m),...
                            'LineWidth',mlineWidth(m),...
                            'Marker',mlineMarker{m});
    end
    daspect([1,1,1]);
    title('Vector Search Space')
    subplot(2,2,2);hold('on')
end



%% Double Check which direction the head is pointed relative to
%  the direction of movement
if dispFigs,
    hfig = figure(230302030);
    set(hfig,'units','centimeters')
    set(hfig,'Position',[11 0 32 25])
    set(hfig,'PaperPositionMode','auto');

    sfac = factor(numel(sessionIds));

    while numel(sfac)>2,
        sfac = sort(sfac,'descend');
        sfac(end-1) = sfac(end)*sfac(end-1);
        sfac(end) = [];
    end
    if numel(sfac)==1
        sfac = [sfac,1];
    end

    
    
    for ses = 1:numel(sessionIds),

        ind = data(:,11)==sessionIds(ses);
        cxyz = cat(2,xyz(ind,[1,2,3],:),permute(data(ind,[10,8,9])+data(ind,[4,2,3]),[1,3,2]));
        %cxyz = cat(2,xyz(ind,[1,2],:),bvecs(:,1,:));

        nframe = size(cxyz,1); %nframe: number of frames (time)
        nmar   = size(cxyz,2);  %nmar: number markers 
        ndim   = size(cxyz,3); %ndim: number of spatial dimensions (xyz)
        j =1:nmar;
        diffMat = permute(cat(4,permute(reshape(repmat(cxyz(:,:,1)',nmar,1)-cxyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(cxyz(:,:,2)',nmar,1)-cxyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(cxyz(:,:,3)',nmar,1)-cxyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);


    % Create Angles between all markers
    % ang(time,marker1,marker2,dim) 
    %     dim: spherical coordinates (theta, phi   , rho     )
    %                                (yaw  , pitch , distance)
        ang = zeros(size(cxyz,1),size(cxyz,2),size(cxyz,2),3);
        for i=1:size(cxyz,2),
            for j=1:size(cxyz,2),
                if i==j,continue,end                    
                switch size(cxyz,3)
                  case 3
                    tang =cell(1,3);
                    [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
                  case 2
                    tang =cell(1,2);
                    [tang{:}] = cart2pol(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
                end
                ang(:,i,j,:) = cell2mat(tang);
            end
        end
        ang(ang(:,1,2,2)~=0,1,1,1)=1;


        vind = nniz(vxy(ind,2));
        tvxy = vxy(ind,2);

        subplot(sfac(2),sfac(1),ses)
        hist2([circ_dist(ang(vind,3,2,1),ang(vind,1,4,1)),tvxy(vind)],40,linspace(-.5,2,40));
        title(['session: ',num2str(ses)])
    end
    print(hfig,'-depsc2',fullfile(odir,ofilename,[ofilename,'summary.eps']));
    print(hfig,'-dpng',fullfile(odir,ofilename,[ofilename,'summary.png']));
end 





% $$$ cxyz = cat(2,xyz(ind,[1],:),bvecs,xyz(ind,[6,7,8],:));
% $$$ c = [0,1,0;...
% $$$      0.5,0,0.5;...
% $$$      1,0,0;...
% $$$      0,0,1;...
% $$$      0,0,1;...
% $$$      0,0,1];
% $$$ 
% $$$ 
% $$$ 
% $$$ cind = 1000;
% $$$ figure
% $$$ hax = axes;
% $$$ marker_options = {};
% $$$ marker_options.style ='o' ;
% $$$ marker_options.size = 8 ;
% $$$ marker_options.erase = 'none'; %{none|xor|background}
% $$$ markers = repmat({0},1,size(cxyz,2));
% $$$ 
% $$$ for l=1:size(cxyz,2),
% $$$     for b=l:size(cxyz,2),
% $$$         sticks{l,b} =animatedline([cxyz(cind(end),l,1),...
% $$$                                     cxyz(cind(end),b,1)],...
% $$$                                 [cxyz(cind(end),l,2),...
% $$$                                     cxyz(cind(end),b,2)],...
% $$$                                 [cxyz(cind(end),l,3),...
% $$$                                     cxyz(cind(end),b,3)]);
% $$$        set(sticks{l,b},'Parent',hax,...
% $$$                      'LineWidth', 2, ... stick_options.width,
% $$$                      'Visible','on',...
% $$$                      'Parent',hax);
% $$$         
% $$$     end
% $$$ end
% $$$ 
% $$$ % $$$ c = [0,1,0;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      1,0,0;...
% $$$ % $$$      0.5,0,0.5;...
% $$$ % $$$      0,1,0;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      1,0,0];
% $$$ 
% $$$ for l=1:size(cxyz,2),
% $$$     markers{l} =animatedline([cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
% $$$                              [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
% $$$                              [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
% $$$     set(markers{l},'Parent',hax, ...
% $$$                    'Marker'         , marker_options.style, ...
% $$$                    'MarkerEdgeColor', c(l,:), ...
% $$$                    'MarkerSize'     , marker_options.size, ...
% $$$                    'MarkerFaceColor', c(l,:), ...
% $$$                    'Visible','on');
% $$$     
% $$$ end    
% $$$ 
% $$$ 
% $$$ xlim([-0.6, 0.6])
% $$$ ylim([-0.5, 0.5])
% $$$ zlim([0, 0.6])
% $$$ daspect([1,1,1])
% $$$ 
% $$$ for cind = 1:10:10000
% $$$     for l=1:size(cxyz,2),
% $$$         for b=l:size(cxyz,2),
% $$$ 
% $$$             clearpoints(sticks{l,b});
% $$$             addpoints(sticks{l,b},[cxyz(cind(end),l,1),...
% $$$                                    cxyz(cind(end),b,1)],...
% $$$                                   [cxyz(cind(end),l,2),...
% $$$                                    cxyz(cind(end),b,2)],...
% $$$                                   [cxyz(cind(end),l,3),...
% $$$                                    cxyz(cind(end),b,3)]);
% $$$             
% $$$         end
% $$$         clearpoints(markers{l});
% $$$         addpoints(markers{l},[cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
% $$$                              [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
% $$$                              [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
% $$$     
% $$$     end
% $$$     drawnow
% $$$     pause(0.001);
% $$$ end


% $$$ 
% $$$ 
% $$$ dispFigs = true;
% $$$ 
% $$$ if dispFigs,
% $$$     set(0,'defaultAxesFontSize',10,...
% $$$           'defaultTextFontSize',10)
% $$$     mkdir(fullfile(odir,ofilename));
% $$$ end
% $$$ 
% $$$ 
% $$$ xyzSampleRate = 12;
% $$$ 
% $$$ 
% $$$ % euler ang vector
% $$$ eang = deg2rad(data(:,[5,6,7]));
% $$$ 
% $$$ 
% $$$ % RxRyRz
% $$$ rMat = cell([size(eang,1),1]);
% $$$ for i = 1:size(eang,1),
% $$$ rMat{i} = [cos(pitch(i)).*cos(yaw(i)),...
% $$$           -cos(pitch(i)).*sin(yaw(i)),...
% $$$            sin(pitch(i));...
% $$$                ...
% $$$            cos(roll(i)).*sin(yaw(i))+sin(roll(i)).*sin(pitch(i)).*cos(yaw(i)),...
% $$$            cos(roll(i)).*cos(yaw(i))-sin(roll(i)).*sin(pitch(i)).*sin(yaw(i)),...
% $$$           -sin(roll(i)).*cos(pitch(i));...
% $$$                ...
% $$$           sin(roll(i)).*sin(yaw(i))-cos(roll(i)).*sin(pitch(i)).*cos(yaw(i)),...
% $$$           sin(roll(i)).*cos(yaw(i))+cos(roll(i)).*sin(pitch(i)).*sin(yaw(i)),...
% $$$           cos(roll(i)).*cos(pitch(i))];
% $$$ end
% $$$ 


