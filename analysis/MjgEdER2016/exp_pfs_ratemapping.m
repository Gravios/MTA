
Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial.load('stc','NN0317R');

OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = 'exp_pfs_ratemapping';
mkdir(fullfile(OwnDir,FigDir))

pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;

%unitt = select_units(Trial,18);
%units = units(Trial.nq.SNR(units)>.75);
Trial.load('nq');
units = find(Trial.nq.AmpSym<0&mrt>1)';


states = {'rear','loc','lloc','hloc','pause','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
nsts = numel(states);


for s = 1:nsts,        
    % Load/Create place fields
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{s};
    defargs = struct2varargin(defargs);        
    pfs{s} = MTAAknnpfs_bs(Trial,defargs{:});      
end


rmcf = [];
for unit = units,
    rm = reshape(cat(3,pfs{3}.plot(unit,'mean'),pfs{4}.plot(unit,'mean')),[],2);
    ind = nniz(rm);
    if sum(ind)>10,
        cf = corr(rm(ind,:));
        rmcf(unit==units) = cf(2);
    else
        rmcf(unit==units) = 0;
    end
end

%figure,
hold on
plot(sort(rmcf))