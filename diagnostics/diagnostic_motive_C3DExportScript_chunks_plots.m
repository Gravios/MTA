itf = c3dserver;

c3dlist = {'E:\motive\ER11-20150520-cof3d_0001-5000.c3d',...
           'E:\motive\ER11-20150520-cof3d_0001-5000_vicon.c3d',...
           'E:\motive\ER11-20150520-cof3d_5001-10000.c3d',...
           'E:\motive\ER11-20150520-cof3d_5001-10000_vicon.c3d'...
           };

       
figure,hold on
shift = 0.001;
for c = 1:numel(c3dlist)
openc3d(itf,0,c3dlist{c});
[xyz,~] = getAll3dTargets(itf);
rescale = 1;
t = regexp(c3dlist{c},'-cof3d_(\d{1,5})-\d{1,5}.+\.','tokens');
t = linspace([cellfun(@str2num,t{1}),cellfun(@str2num,t{1})+size(xyz,1),size(xyz,1)]);
if ~isempty(regexp(c3dlist{c},'_vicon'))
    rescale = 1000;
end
plot(t,xyz(:,1,1)./rescale+shift*c)
closec3d(itf)
end

legend(cellfun(@(c)c{1},regexp(c3dlist,'cof3d_\d{1,5}-\d{1,5}.+\.','match'), 'UniformOutput', false));


itf.release