function PHI=phase(G,varargin)
%PHASE  Computes the phase of a complex vector
%
%   PHI=phase(G)
%
%   G is a complex-valued row vector and PHI is returned as its
%   phase (in radians), with an effort made to keep it continuous
%   over the pi-borders.

%   L. Ljung 10-2-86
%   Copyright 1986-2004 The MathWorks, Inc.
%   $Revision: 1.5.4.2 $  $Date: 2004/07/31 23:24:49 $
[flag] = DefaultArgs(varargin,{1,''});

ndim = ones(1,4);
ndim(1:size(size(G)))=size(G);
G = G(:);
PHI = atan2(imag(G),real(G));
PHI = reshape(PHI,ndim);

if strcmp('unwrap',flag),
    for j = 1:ndim(2),
        for k = 1:ndim(3),
            DF=PHI(1:ndim(1)-1,j,k)-PHI(2:ndim(1),j,k);
            I=find(abs(DF)>3.5);
            if ~isempty(I)
                for i=I
                    PHI(:,j,k) = PHI(:,j,k)+2*pi*sign(DF(i))*[zeros(1,i),ones(1,ndim(1)-i)];
                end
            end
        end
    end
end

end