
% req20171215 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: use placefield rate and area to search for point 
%               spatial representation
%  Bugs: NA


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
Trial = MTATrial.validate('Ed10-20140817.cof.gnd');

state = 'pause&theta';

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
pargs.overwrite = true;
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





i = [-100:10:100];
j = [-100:10:100];
k = [-50:10:50];

[mrt,mrx] = pft{21,21}.maxRate(units);

ru = reshape(unity(reshape(r,[],numel(units))),size(r));
ind = sqrt(sum(mrx.^2,2))<350;
ruo = MTADfet.encapsulate(Trial,mean(ru(:,:,:,ind),4),1,'pmin','pmin','p');
ruo.filter('RectFilter',5,8);
ruo.data = permute(ruo.data,[2,3,1]);
ruo.filter('RectFilter',5,8);
ruo.data = permute(ruo.data,[2,3,1]);
ruo.filter('RectFilter',3,8);
ruo.data = permute(ruo.data,[2,3,1]);

[ruomin,ruoval] = LocalMinimaN(-ruo.data,0,[3,3,3]);

hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hax = gobjects([1,numel(k)]);
ruolmin = nan([numel(k),2]);
for z = fliplr(1:numel(k)),
    hax(z) = subplot(numel(k),1,numel(k)-z+1);
    hold('on');
    imagesc(i,j,-ruo(:,:,z)');
    ruolmin(z,:) = LocalMinima2(-ruo(:,:,z),10,15);
    scatter(i(ruolmin(z,1)),j(ruolmin(z,2)),20,'m','filled');
    if any(ruomin(:,3)==z)
        scatter(i(ruomin(ruomin(:,3)==z,1)),j(ruomin(ruomin(:,3)==z,2)),40,'w');
    end
    title([num2str(k(z)),' mm']);
    xlim([-100,100]);
    ylim([-100,100]);
end
hfig.Position = [0.5,0.5,6,40];
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
FigName = ['pfs_optimal_pos_head_3d_',Trial.filebase];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));


hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hax = gobjects([1,2]);
hax(1) = subplot(121);
plot(i(ruolmin(:,1)),k,'.');
xlabel('longitudinal shift mm')
ylabel('vertical shift mm')
hax(2) = subplot(122);
plot(j(ruolmin(:,2)),k,'.')
xlabel('lateral shift mm')
ylabel('vertical shift mm')
hfig.Position = [0.5,0.5,10,6];
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
FigName = ['pfs_optimal_pos_head_3d_maxima',Trial.filebase];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));



% $$$ 
% $$$ 
% $$$ FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ offset = 0;
% $$$ hax = gobjects([12,3]);
% $$$ for u = 1:12,
% $$$ % PLOT placefields
% $$$     hax(u,1) = subplot2(12,18,u,[1:4]);
% $$$     plot(pft{11,6,6},units(u+offset));
% $$$     if u~=12, clear_axes_ticks(hax(u,1));
% $$$     else,     xlabel('mm');
% $$$               ylabel('mm');
% $$$     end
% $$$     if u==1, title('placefield'), end
% $$$ % PLOT max rate as function of head translation
% $$$     hax(u,2) = subplot2(12,18,u,[5:11]);
% $$$     imagesc(i,k,sq(r(:,6,:,u+offset))');axis('xy');
% $$$     if u~=12, clear_axes_ticks(hax(u,2));  
% $$$     else,     hax(u,2).YTickLabelMode = 'manual';
% $$$               hax(u,2).YTickLabel = {};
% $$$               hax(u,2).XTickMode = 'manual';
% $$$               hax(u,2).XTick = [-100,-50,0,50,100];
% $$$               xlabel('mm');
% $$$     end
% $$$     if u==1, title('max rate'), end
% $$$ % PLOT half-max-rate area as function of head translation
% $$$     hax(u,3) = subplot2(12,18,u,[12:18]);
% $$$     imagesc(i,k,sq(p(:,6,:,u+offset))');axis('xy');
% $$$     if u~=12, clear_axes_ticks(hax(u,3));
% $$$     else,     hax(u,3).YAxisLocation = 'right';
% $$$               hax(u,3).XTickMode = 'manual';
% $$$               hax(u,3).XTick = [-100,-50,0,50,100];
% $$$               xlabel('mm');
% $$$               ylabel('mm');
% $$$     end
% $$$     if u==1, title('half-max-rate area'); end
% $$$ end
% $$$ suptitle('max rate and half-max-rate area as function of head translation');
% $$$ hfig.Position = [0.5,0.5,16,40];
% $$$ af(@(h)  set(h,'Units','centimeters'),  hax);
% $$$ af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
% $$$ af(@(h)  set(h,'Position',[h.Position(1:2),4,2]),  hax(:,2:3));
% $$$ FigName = ['pfs_rate_max_area_min','_head_xz_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
% $$$ print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));

% $$$ FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ offset = 0;
% $$$ hax = gobjects([12,3]);
% $$$ for u = 1:12,
% $$$ % PLOT placefields
% $$$     hax(u,1) = subplot2(12,18,u,[1:4]);
% $$$     plot(pft{11,6,6},units(u+offset));
% $$$     if u~=12, clear_axes_ticks(hax(u,1));
% $$$     else,     xlabel('mm');
% $$$               ylabel('mm');
% $$$     end
% $$$     if u==1, title('placefield'), end
% $$$ % PLOT max rate as function of head translation
% $$$     hax(u,2) = subplot2(12,18,u,[5:11]);
% $$$     imagesc(i,j,sq(r(:,:,6,u+offset))');axis('xy');
% $$$     if u~=12, clear_axes_ticks(hax(u,2));  
% $$$     else,     hax(u,2).YTickLabelMode = 'manual';
% $$$               hax(u,2).YTickLabel = {};
% $$$               hax(u,2).XTickMode = 'manual';
% $$$               hax(u,2).XTick = [-100,-50,0,50,100];
% $$$               xlabel('mm');
% $$$     end
% $$$     if u==1, title('max rate'), end
% $$$ % PLOT half-max-rate area as function of head translation
% $$$     hax(u,3) = subplot2(12,18,u,[12:18]);
% $$$     imagesc(i,j,sq(p(:,:,6,u+offset))');axis('xy');
% $$$     if u~=12, clear_axes_ticks(hax(u,3));
% $$$     else,     hax(u,3).YAxisLocation = 'right';
% $$$               hax(u,3).XTickMode = 'manual';
% $$$               hax(u,3).XTick = [-100,-50,0,50,100];
% $$$               xlabel('mm');
% $$$               ylabel('mm');
% $$$     end
% $$$     if u==1, title('half-max-rate area'); end
% $$$ end
% $$$ suptitle('max rate and half-max-rate area as function of head translation');
% $$$ hfig.Position = [0.5,0.5,16,40];
% $$$ af(@(h)  set(h,'Units','centimeters'),  hax);
% $$$ af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
% $$$ af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:3));
% $$$ FigName = ['pfs_rate_max_area_min','_head_xy_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
% $$$ print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));







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
pargs.states = state;
pargs.numIter = 1;
pargs.overwrite = true;
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
offset = 12;
hax = gobjects([12,3]);
clf();
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



% search space
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


ny = 12; 
for offset = 0:ny:numel(units)-1,
    clf();
    hax = gobjects([1,4]);
    for u = 1:ny,
        if numel(units) < u+offset,break;end;        

% PLOT accg
        h=1;
        hax(u,h) = subplot2(ny,18,u,[1:2]);
        bar(tbins,accg(:,units(u+offset)));
        xlim(tbins([1,end]));        
        if u~=ny, clear_axes_ticks(hax(u,h));  end
        title(num2str(units(u+offset)));

% PLOT placefield
        h = h+1;
        hax(u,h) = subplot2(ny,18,u,[3:4]);
        plot(pfd{21,21},units(u+offset));
        if u~=ny, clear_axes_ticks(hax(u,h));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('placefield'), end

        
% PLOT DRZ field pitch 
        h = h+1;
        if numel(units) < u+offset,break;end;
        hax(u,h) = subplot2(ny,18,u,[5:6]);
        plot(dfs{1},units(u+offset),'mazeMaskFlag',false);
        if u~=ny, clear_axes_ticks(hax(u,h));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('drzfield pitch'), end

% PLOT DRZ field height         
        h = h+1;        
        if numel(units) < u+offset,break;end;
        hax(u,h) = subplot2(ny,18,u,[7:8]);
        plot(dfs{2},units(u+offset),'mazeMaskFlag',false);
        if u~=ny, clear_axes_ticks(hax(u,h));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('drzfield Height'), end        

% PLOT DRZ field height         
        h = h+1;        
        if numel(units) < u+offset,break;end;
        hax(u,h) = subplot2(ny,18,u,[9:10]);
        plot(dfs{3},units(u+offset),'mazeMaskFlag',false);
        if u~=ny, clear_axes_ticks(hax(u,h));
        else,     xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('drzfield RHM'), end        

        
        
        
        % PLOT max rate as function of 3d head translation XZ
        h = h+1;
        hax(u,h) = subplot2(ny,18,u,[11:12]);
        imagesc(id,kd,sq(r(:,round(numel(jd)/2),:,u+offset))');axis('xy');
        if u~=ny, clear_axes_ticks(hax(u,h));  
        else,     hax(u,h).YTickLabelMode = 'manual';
            hax(u,h).YTickLabel = {};
            hax(u,h).XTickMode = 'manual';
            hax(u,h).XTick = [-100,-50,0,50,100];
            xlabel('mm');
        end
        if u==1, title('3d XZ-plane'), end        

% PLOT max rate as function of 3d head translation XY
        h = h+1;
        hax(u,h) = subplot2(ny,18,u,[13:14]);
        imagesc(id,jd,sq(r(:,:,round(numel(kd)/2),u+offset))');axis('xy');
        if u~=ny, clear_axes_ticks(hax(u,h));  
        else,     hax(u,h).YTickLabelMode = 'manual';
            hax(u,h).YTickLabel = {};
            hax(u,h).XTickMode = 'manual';
            hax(u,h).XTick = [-100,-50,0,50,100];
            xlabel('mm');
        end    
        if u==1, title('3d XY-plane'), end    
        
% PLOT max rate as function of 2d head translation XY
        h = h+1;
        hax(u,h) = subplot2(ny,18,u,[15:16]);
        imagesc(i,j,sq(rd(:,:,u+offset))');axis('xy');
        if u~=ny, clear_axes_ticks(hax(u,h));  
        else,     hax(u,3).YAxisLocation = 'right';
            hax(u,3).XTickMode = 'manual';
            hax(u,3).XTick = [-100,-50,0,50,100];
            xlabel('mm');
            ylabel('mm');
        end
        if u==1, title('2d XY-plane'), end        
    end
    suptitle(['max rate as function of head translation: ',Trial.filebase]);
    hfig.Position = [0.5,0.5,28,40];
    af(@(h)  set(h,'Units','centimeters'),  hax);
    af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
    af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:end));
    FigName = ['pfs_rate_max_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
    print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));
end




figure,plot(i,mean(bsxfun(@rdivide,rd(:,21,:),mean(rd(:,21,:))),3))

[mrt,mrx] = pft{21,21}.maxRate(units);


rd = reshape(unity(reshape(rd,[],numel(units))),size(rd));

ind = sqrt(sum(mrx.^2,2))<250;
%figure,imagesc(mean(rd(:,:,ind),3)');

rdo = MTADfet.encapsulate(Trial,mean(rd(:,:,ind),3),1,'pmin','pmin','p');
rdo.filter('RectFilter',5,5);
rdo.data = rdo.data';
rdo.filter('RectFilter',5,5);
rdo.data = rdo.data';



rdmin = LocalMinima2(-rdo.data,0,15);
figure,hold('on');
imagesc(i,j,-rdo.data');
scatter(i(rdmin(1)),j(rdmin(2)),20,'m','filled');
axis('xy');
axis('tight');
ylabel('mm');
xlabel('mm');
title('Placefield optimal marker position relative to head');
FigName = ['pfs_optimal_pos_head_',Trial.filebase];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));




rumin = LocalMinimaN(-ruo.data,0,[10,10,5]);
figure,hold('on');
imagesc(id,jd,-ruo.data(:,:,6)');
scatter(id(rumin(1)),jd(rumin(2)),20,'m','filled');
axis('xy');
axis('tight');
ylabel('mm');
xlabel('mm');
title('Placefield optimal marker position relative to head');
FigName = ['pfs_optimal_pos_3d_head_',Trial.filebase];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));
