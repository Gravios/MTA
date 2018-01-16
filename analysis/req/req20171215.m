
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
Trial = MTATrial.validate('jg05-20120311.cof.all');

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
j = [-100:10:100];
k = [-50:10:50];

txyz = xyz(:,nhm,:);
ind = nniz(txyz);


pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units = units;
pargs.states = 'theta-groom-sit';
pargs.numIter = 1;
pargs.overwrite = false;
pargs.halfsample = false;
pft = cell([numel(i),numel(j),numel(k)]);
for x = 1:numel(i),
    for y = 1:numel(j),
        for z = 1:numel(k)
            sxyz = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*k(z),txyz);
            pargs.tag = ['ms-x',num2str(i(x)),'y',num2str(j(y)),'z',num2str(k(z))];
            pargs.xyzp = MTADfet.encapsulate(Trial,sxyz(:,1,[1,2]),xyz.sampleRate,'','','');
            pfsArgs = struct2varargin(pargs);                    
            pft{x,y,z} = MTAApfs(Trial,pfsArgs{:});
        end
    end
end

mrates = cf(@(p,u)  p.maxRate(u),  pft,repmat({units},size(pft)));
pareas = cf(@(p,m)  permute(sum(bsxfun(@gt,p.data.rateMap,m'./2)),[1,3,4,2]),  pft,  mrates);
p = cell2mat(pareas);

r = cf(@(m) permute(m,[2,3,4,1]),  mrates);
r = cell2mat(r);




FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
offset = 0;
hax = gobjects([12,3]);
for u = 1:12,
% PLOT placefields
    hax(u,1) = subplot2(12,18,u,[1:4]);
    plot(pft{11,6,6},units(u+offset));
    if u~=12, clear_axes_ticks(hax(u,1));
    else,     xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('placefield'), end
% PLOT max rate as function of head translation
    hax(u,2) = subplot2(12,18,u,[5:11]);
    imagesc(i,k,sq(r(:,6,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,2));  
    else,     hax(u,2).YTickLabelMode = 'manual';
              hax(u,2).YTickLabel = {};
              hax(u,2).XTickMode = 'manual';
              hax(u,2).XTick = [-100,-50,0,50,100];
              xlabel('mm');
    end
    if u==1, title('max rate'), end
% PLOT half-max-rate area as function of head translation
    hax(u,3) = subplot2(12,18,u,[12:18]);
    imagesc(i,k,sq(p(:,6,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,3));
    else,     hax(u,3).YAxisLocation = 'right';
              hax(u,3).XTickMode = 'manual';
              hax(u,3).XTick = [-100,-50,0,50,100];
              xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('half-max-rate area'); end
end
suptitle('max rate and half-max-rate area as function of head translation');
hfig.Position = [0.5,0.5,16,40];
af(@(h)  set(h,'Units','centimeters'),  hax);
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
af(@(h)  set(h,'Position',[h.Position(1:2),4,2]),  hax(:,2:3));
FigName = ['pfs_rate_max_area_min','_head_xz_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));

FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
offset = 0;
hax = gobjects([12,3]);
for u = 1:12,
% PLOT placefields
    hax(u,1) = subplot2(12,18,u,[1:4]);
    plot(pft{11,6,6},units(u+offset));
    if u~=12, clear_axes_ticks(hax(u,1));
    else,     xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('placefield'), end
% PLOT max rate as function of head translation
    hax(u,2) = subplot2(12,18,u,[5:11]);
    imagesc(i,j,sq(r(:,:,6,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,2));  
    else,     hax(u,2).YTickLabelMode = 'manual';
              hax(u,2).YTickLabel = {};
              hax(u,2).XTickMode = 'manual';
              hax(u,2).XTick = [-100,-50,0,50,100];
              xlabel('mm');
    end
    if u==1, title('max rate'), end
% PLOT half-max-rate area as function of head translation
    hax(u,3) = subplot2(12,18,u,[12:18]);
    imagesc(i,j,sq(p(:,:,6,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,3));
    else,     hax(u,3).YAxisLocation = 'right';
              hax(u,3).XTickMode = 'manual';
              hax(u,3).XTick = [-100,-50,0,50,100];
              xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('half-max-rate area'); end
end
suptitle('max rate and half-max-rate area as function of head translation');
hfig.Position = [0.5,0.5,16,40];
af(@(h)  set(h,'Units','centimeters'),  hax);
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:3));
FigName = ['pfs_rate_max_area_min','_head_xy_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));







% 2d placefield optimizaiton search

xyz = Trial.load('xyz','trb');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom', [0.5,1,0.5],{{'head_back','head_front',[0,0,1]}},hcom);
hcom(:,:,3) = 0;
xyz.data(:,:,3) = 0;

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


i = [-100:5:100];
j = [-100:5:100];

txyz = xyz(:,nhm,:);

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units = units;
pargs.states = 'theta-groom-sit';
pargs.numIter = 1;
pargs.overwrite = false;
pargs.halfsample = false;
pfd = cell([numel(i),numel(j)]);
for x = 1:numel(i),
    for y = 1:numel(j),
            sxyz = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*0,txyz);
            pargs.tag = ['ms-x',num2str(i(x)),'y',num2str(j(y))];
            pargs.xyzp = MTADfet.encapsulate(Trial,sxyz(:,1,[1,2]),xyz.sampleRate,'','','');
            pfsArgs = struct2varargin(pargs);
            pfd{x,y} = MTAApfs(Trial,pfsArgs{:});
    end
end

mratesd = cf(@(p,u)  p.maxRate(u),  pfd,repmat({units},size(pfd)));
pareasd = cf(@(p,m)  permute(sum(bsxfun(@gt,p.data.rateMap,m'./2)),[1,3,2]),  pfd,  mratesd);
pd = cell2mat(pareasd);

rd = cf(@(m) permute(m,[2,3,1]),  mratesd);
rd = cell2mat(rd);


FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
offset = 0;
hax = gobjects([12,3]);
for u = 1:12,
% PLOT placefields
    hax(u,1) = subplot2(12,18,u,[1:4]);
    plot(pfd{21,21},units(u+offset));
    if u~=12, clear_axes_ticks(hax(u,1));
    else,     xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('placefield'), end
% PLOT max rate as function of head translation
    hax(u,2) = subplot2(12,18,u,[5:11]);
    imagesc(i,j,sq(rd(:,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,2));  
    else,     hax(u,2).YTickLabelMode = 'manual';
              hax(u,2).YTickLabel = {};
              hax(u,2).XTickMode = 'manual';
              hax(u,2).XTick = [-100,-50,0,50,100];
              xlabel('mm');
    end
    if u==1, title('max rate'), end
% PLOT half-max-rate area as function of head translation
    hax(u,3) = subplot2(12,18,u,[12:18]);
    imagesc(i,j,sq(pd(:,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,3));
    else,     hax(u,3).YAxisLocation = 'right';
              hax(u,3).XTickMode = 'manual';
              hax(u,3).XTick = [-100,-50,0,50,100];
              xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('half-max-rate area'); end
end
suptitle('max rate and half-max-rate area as function of head translation');
hfig.Position = [0.5,0.5,16,40];
af(@(h)  set(h,'Units','centimeters'),  hax);
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:3));
FigName = ['pfs_rate_max_area_min','_2d_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));


id = [-100:10:100];
jd = [-100:10:100];
kd = [-50:10:50];




% LOAD DRZ fields
dfs = cell([1,3]);
[dfs{:}] = MjgER2016_drzfields(Trial,units,false);%,pfstats);
dfst = {'pitch','height','rhm'};

% LOAD unit auto-correlogram
[accg,tbins] = autoccg(MTASession.validate(Trial.filebase));



FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';

for offset = 0:12:numel(units)-1,
    clf();
    hax = gobjects([1,4]);
    for u = 1:12,
        if numel(units) < u+offset,break;end;        

% PLOT accg
        hax(u,1) = subplot2(12,18,u,[1:2]);
        plot(pfd{21,21},units(u+offset));

% PLOT placefield
        hax(u,1) = subplot2(12,18,u,[1:2]);
        plot(pfd{21,21},units(u+offset));
        if u~=12, clear_axes_ticks(hax(u,1));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('placefield'), end

% PLOT DRZ field pitch 
        if numel(units) < u+offset,break;end;
        hax(u,1) = subplot2(12,18,u,[3:4]);
        plot(pfd{1},units(u+offset));
        if u~=12, clear_axes_ticks(hax(u,1));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('placefield'), end

        
        if numel(units) < u+offset,break;end;
        hax(u,1) = subplot2(12,18,u,[3:4]);
        plot(pfd{2},units(u+offset));
        if u~=12, clear_axes_ticks(hax(u,1));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('placefield'), end
        
        
        
        % PLOT max rate as function of 3d head translation XZ
        h = 2;
        hax(u,h) = subplot2(12,18,u,[7:10]);
        imagesc(id,kd,sq(r(:,round(numel(jd)/2),:,u+offset))');axis('xy');
        if u~=12, clear_axes_ticks(hax(u,h));  
        else,     hax(u,h).YTickLabelMode = 'manual';
            hax(u,h).YTickLabel = {};
            hax(u,h).XTickMode = 'manual';
            hax(u,h).XTick = [-100,-50,0,50,100];
            xlabel('mm');
        end
        if u==1, title('3d XZ-plane'), end        
        % PLOT max rate as function of 3d head translation XY
        h = 3;
        hax(u,h) = subplot2(12,18,u,[11:14]);
        imagesc(id,jd,sq(r(:,:,round(numel(kd)/2),u+offset))');axis('xy');
        if u~=12, clear_axes_ticks(hax(u,h));  
        else,     hax(u,h).YTickLabelMode = 'manual';
            hax(u,h).YTickLabel = {};
            hax(u,h).XTickMode = 'manual';
            hax(u,h).XTick = [-100,-50,0,50,100];
            xlabel('mm');
        end    
        if u==1, title('3d XY-plane'), end    
        % PLOT max rate as function of 2d head translation XY
        h = 4;
        hax(u,h) = subplot2(12,18,u,[15:18]);
        imagesc(i,j,sq(rd(:,:,u+offset))');axis('xy');
        if u~=12, clear_axes_ticks(hax(u,h));  
        else,     hax(u,3).YAxisLocation = 'right';
            hax(u,3).XTickMode = 'manual';
            hax(u,3).XTick = [-100,-50,0,50,100];
            xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('2d XY-plane'), end        
    end
    suptitle('max rate as function of head translation');
    hfig.Position = [0.5,0.5,26,40];
    af(@(h)  set(h,'Units','centimeters'),  hax);
    af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
    af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:end));
    FigName = ['pfs_rate_max_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
    print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));
end