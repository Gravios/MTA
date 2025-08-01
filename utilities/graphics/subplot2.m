%function h = subplot2b(ny,nx,py,px)
% 2d indexing of subplot -nx, ny as before
% but px and py are subplot index coordinates 
% e.g. subplot(2,2,3) = subplot2(2,2,2,1);
%
% Mod history 
% JG - added feature: now suports retangular subplots
function h = subplot2(ny,nx,py,px)

sc = 0:ny-1;
p = [];
for i = 1:length(py),
p = cat(2,p,px+sc(py(i)).*nx);
end

h = subplot(ny,nx,p(:),'replace');

