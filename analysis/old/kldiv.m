function Dkl = kldiv(varargin)

%Testing vars
if nargn==0,
    Trial = MTATrial('jg05-20120317');    
    Data = create(MTADang,Trial,Trial.load('xyz').filter('ButFilter',3,3,'low'));
    Data.data = [[diff(Data(:,1,4,2));0].*Data.sampleRate,Data(:,3,4,2)];
    varargin{1} = Data;
    varargin{2} = {'gper-rear','rear'};
    varargin{3} = {[-4,4],[-1,pi/2]};

    xyz = Trial.load('xyz').filter('ButFilter',3,3,'low')
    ang =  create(MTADang,Trial,xyz);
    Data = Trial.xyz.copy;
    Data.data = [ang(:,1,4,2)];
    varargin{1} = Data;
    varargin{2} = {'gper-rear','rear'};
    varargin{3} = {[.3,pi/2]};
    varargin{4} = [1500];
end
    
defargs = {...  Data
               [],...
           ...  States
               {'gper-rear','rear'},...
           ...  lims 
               {[-2,2],[-2,2]},...
           ...  nbins
               [100,100]};

[Data,states,lims,nbins] = DefaultArgs(varargin,defargs,true);

ndims = Data.size(2);

assert(numel(lims)==numel(nbins)&&numel(lims)==ndims,['MTA:' ...
                    'analysis:kldiv second dimension Data must equal ' ...
                    'the number of elements in bins and lims'])


lims = reshape(cell2mat(lims),2,[])';
Bins = cell(1,ndims);


Occupancy = {};
for s = states,
    pos = Data(Trial.stc{s{:}},:);
    %pos = Data(s,:);    

    k = nbins'./abs(diff(lims,1,2));
    msize = [round(abs(diff(lims,1,2)).*k);1];
    for i = 1:ndims
        Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
    end

    Pos = round((pos-repmat(lims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
    for i = 1:ndims   
        Pos(Pos(:,i)<1|Pos(:,i)>nbins(i)|~nniz(Pos),:) = [];
    end

    Occupancy{end+1} = accumarray(Pos,1,msize');
    Occupancy{end} = Occupancy{end}/sum(Occupancy{end}(:));

end

Dkl = Occupancy{1}(:).*log(Occupancy{1}(:)./Occupancy{2}(:));
Dkl = sum(Dkl(nniz(Dkl)));

ldkl(end+1) = Dkl;

end