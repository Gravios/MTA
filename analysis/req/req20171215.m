
% $$$ % SEARCH for position which minimizes placefield size
% $$$ 
% $$$ units = Trial.spk.map(:,1)';
% $$$ markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'};
% $$$ for m = 1:numel(markers),
% $$$     Trial.trackingMarker = markers{m};
% $$$     defargs = get_default_args('MjgER2016','MTAApfs','struct');
% $$$     defargs.units = units;
% $$$     defargs.overwrite = true;
% $$$     defargs.states = 'theta-groom-sit';
% $$$     defargs = struct2varargin(defargs);        
% $$$     pft{m} = MTAApfs(Trial,defargs{:});      
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ % place fields generated from each marker along the spine and head
% $$$ figure();
% $$$ for unit = units,
% $$$     clf();
% $$$     for m = 1:numel(markers),
% $$$         subplot(1,numel(markers),m);
% $$$         plot(pft{m},unit);    
% $$$         title([num2str(unit),' ',markers{m}]);
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 
% $$$ % normalized spatial information as function of rostro-caudal distance
% $$$ figure();
% $$$ plot(bsxfun(@rdivide,cat(2,e{:})',max(cat(2,e{:})')));

Trial = MTATrial.validate('jg05-20120310.cof.all');

units = select_placefields(Trial);

% Now generate a head local gradient of spatial information around the head

xyz = Trial.load('xyz','trb');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom', [0.5,1,0.5],{{'head_back','head_front',[0,0,1]}},hcom);

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
nm = nz.*20+hcom;
xyz.addMarker('htx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;
xyz.addMarker('hbx',  [0.5,1,0.5],[],nm);
nhm = {'hcom'};%,'hbx','hrx','htx'};


i = [-100:10:100];
j = [-50:10:50];
k = [-50:10:50];

txyz = xyz(:,nhm,:);
ind = nniz(txyz);


pft = cell([numel(i),numel(j),numel(k)]);
for x = 1:numel(i),
    for y = 1:numel(j),
        for z = 1:numel(k)
            sxyz = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*k(z),txyz);
            defargs = get_default_args('MjgER2016','MTAApfs','struct');
            defargs.units = units;
            defargs.states = 'theta-groom-sit';
            defargs.numIter = 1;
            %defargs.overwrite = true;
            defargs.halfsample = false;
            defargs.tag = ['ms-x',num2str(i(x)),'y',num2str(j(y)),'z',num2str(k(z))];
            defargs.xyzp = MTADfet.encapsulate(Trial,sxyz(:,1,[1,2]),xyz.sampleRate,'','','');
            defargs = struct2varargin(defargs);        
            
            pft{x,y,z} = MTAApfs(Trial,defargs{:});
        end
    end
end

msi = cell2mat(cf(@(p)  mean(p.data.si,'omitnan'), pft));
ssi = cell2mat(cf(@(p)  std(p.data.si,'omitnan'),  pft));

figure,
for x = 2:2:numel(i),
    subplot2(2,10,1,x/2);
    imagesc(j,k,sq(msi(x,:,:))')
    subplot2(2,10,2,x/2);    
    imagesc(j,k,sq(ssi(x,:,:))')
end

%sxyz(nniz(sxyz),1,:) = ButFilter(sxyz(nniz(sxyz),1,:),3,2/(xyz.sampleRate/2),'low');



% LOAD placefields MTAApfs THETA
defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'theta-groom-sit';
defargs = struct2varargin(defargs);        
pft = MTAApfs(Trial,defargs{:});      

xyz = preproc_xyz(Trial,'LOAD_TRB_XYZ');

ang = create(MTADang,Trial,xyz);


sts = 'rwnpms';
figure,
edx = linspace(30,60,50);
edy = linspace(20,80,50);
for s = 1:numel(sts),
    ind = stc{sts(s)};
    subplot(2,3,s);
    hist2([ang(ind,'spine_middle','spine_upper',3),ang(ind,'spine_upper','hcom',3)],edx,edy);
    grid('on');
end

