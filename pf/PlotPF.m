function [RateMap, Bins, MRate, SI, Spar] = PlotPF(Session,spkpos,pos,varargin)
[binDims,Smooth,type] = DefaultArgs(varargin,{50,[],'xy'});

Nbin = diff(Session.maze.boundaries(1:numel(binDims),:),1,2)./binDims';

%% Constraint to maze is forced
switch type
  case 'xy'
    Xmin = Session.maze.boundaries(1,1);
    Xmax = Session.maze.boundaries(1,2);
    Ymin = Session.maze.boundaries(2,1);
    Ymax = Session.maze.boundaries(2,2);

    %% scaling factor for rounding position
    dx = Xmax - Xmin; 
    dy = Ymax - Ymin; 
    k = [Nbin(1)/dx Nbin(2)/dy];

    %% matrix size
    msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]);

    Bin1 = ([1:msize(1)]-1)/k(1) + Xmin+round(k(1)^-1/2);
    Bin2 = ([1:msize(2)]-1)/k(2) + Ymin+round(k(2)^-1/2);    
    Bins = [Bin1(:),Bin2(:)];

   case 'xyz'
     Xmin = Session.Maze.boundaries(1,1);
     Xmax = Session.Maze.boundaries(1,2);
     Ymin = Session.Maze.boundaries(2,1);
     Ymax = Session.Maze.boundaries(2,2);
     Zmin = Session.Maze.boundaries(2,1);
     Zmax = Session.Maze.boundaries(2,2);
 
     %% scaling factor for rounding position
     dx = Xmax - Xmin; 
     dy = Ymax - Ymin; 
     dz = Zmax - Zmin; 
     k = Nbin./[dx;dy;dz];
 
     %% matrix size
     msize = round(sum(abs([Xmin,Xmax;Ymin,Ymax;Zmin,Zmax]))'.*k);
 
     Bin1 = ([1:msize(1)]-1)/k(1) + Xmin+round(k(1)^-1/2);
     Bin2 = ([1:msize(2)]-1)/k(2) + Ymin+round(k(2)^-1/2);    
     Bin3 = ([1:msize(2)]-1)/k(3) + Ymin+round(k(3)^-1/2);    
     Bins = [Bin1(:),Bin2(:),Bin3(:)];    

  case 'pfcrz'
    switch Session.Maze.shape
      case 'circle'
        Xmin = 0;
        Xmax = Session.Maze.boundaries(1,2)-Session.Maze.boundaries(1,1);
        Ymin = Session.Maze.boundaries(3,1);
        Ymax = Session.Maze.boundaries(3,2);

        %% scaling factor for rounding position
        dx = Xmax - Xmin; 
        dy = Ymax - Ymin; 
        k = [Nbin/dx Nbin/dy];

        %% matrix size
        msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]);

        Bin1 = ([1:msize(1)]-1)/k(1);
        Bin2 = ([1:msize(2)]-1)/k(2);
        Bins = [Bin1(:),Bin2(:)];

      case 'square'
    end

end



%% rounded position
X = round((pos(:,1)-Xmin)*k(1))+1;
Y = round((pos(:,2)-Ymin)*k(2))+1;

%% Push back in any stray bins
X(X>Nbin(1)) = Nbin(1);
Y(Y>Nbin(2)) = Nbin(2);
X(X<1) = 1;
Y(Y<1) = 1;


%% Occupancy
Occ = Accumulate([X Y],1,msize)./Session.xyz.sampleRate;

%% spike count
spikep(:,1) = round((spkpos(:,1)-Xmin)*k(1))+1;
spikep(:,2) = round((spkpos(:,2)-Ymin)*k(2))+1;
Count = Accumulate(spikep,1,msize);

%% smooth
if isempty(Smooth)
  Smooth = Nbin/3000;
end

r1 = (-msize(1):msize(1))/msize(1);
r2 = (-msize(2):msize(2))/msize(2);


% $$$ [x,y,z] = meshgrid(-msize(1):msize(1),-msize(2):msize(2),-msize(3):msize(3));
% $$$ [x,y,z] = meshgrid([-25:25]);
% $$$ Smoother = exp(-x.^2/Smooth^2/2-y.^2/Smooth^2/2-z.^2/Smooth^2/2);
% $$$ Smoother = Smoother./sum(Smoother(:));



Smoother1 = exp(-r1.^2/Smooth^2/2);
Smoother1 = Smoother1/sum(Smoother1);
Smoother2 = exp(-r2.^2/Smooth^2/2);
Smoother2 = Smoother2/sum(Smoother2);

SCount = conv2(Smoother1,Smoother2,Count,'same');
SOcc = conv2(Smoother1,Smoother2,Occ,'same');


OccThresh = 0.06;
%% Find the total occupancy and each pixels 
%% probability of occupancy
gtind = SOcc>OccThresh;
TotOcc = sum(SOcc(gtind));
POcc = SOcc./TotOcc;
%% Rate Map
RateMap = NaN(prod(Nbin),1);
RateMap(gtind) = SCount(gtind)./SOcc(gtind);
RateMap = RateMap(:);
RateMap(~gtind) = NaN;
%% Find the units overall mean rate given the 
%% current state
MRate = sum(SCount(gtind))/TotOcc;
SI = nansum(POcc(gtind).*(RateMap(gtind)./MRate).*log2(RateMap(gtind)./MRate));
Spar = 1/nansum(POcc(gtind).*RateMap(gtind).^2./MRate.^2);

