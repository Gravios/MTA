

%Segmentation JPDF walk
xyz = Trial.load('xyz');
vxy = xyz.copy;
vxy = vxy.vel([1,7],[1,2]);
vxy.filter('ButFilter',3,2.4);
vxy.data(vxy.data<1e-3) = 1e-3;
tag = 'speed_SpineL_vs_HeadF';

hedgs    = {linspace(-1,2,75)};
hedgs(2) = {linspace(-1,2,75)};
edgs    = {linspace(-1,2,75)};
edgs(2) = {linspace(-1,2,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});

sts = 'wrnpms';
stc = [0,0,1;...
       1,0,0;...
       0,1,0;...
       0,1,1;...
       1,0,1;...
       1,1,0];
figure, hold on
hist2(log10(vxy(Trial.stc{'a'},:)),edgs{1},edgs{2});

for i = 1:numel(sts),
    b = log10(vxy(Trial.stc{sts(i)},:));
    o = hist2(b,hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[30,30],'linewidth',2.5,'Color',stc(i,:))
end
legend(states)

xlabel('log10(xy_{speed}) of Lower Body');
ylabel('log10(xy_{speed}) of Head');