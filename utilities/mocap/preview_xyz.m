function hfig = preview_xyz(c3dpath,hfig,hax)

%hfig = [];
if isempty(hfig),
    hfig = figure();
end
if isempty(hax)
    hax = axes();
end

axes(hax);
hold(hax,'on');

files = dir(fullfile(c3dpath,'*.c3d.mat'));
for f = files'
    load(fullfile(c3dpath,f.name));
    plot3(xyzpos(:,2,1),xyzpos(:,2,2),xyzpos(:,2,3),'.');
end





