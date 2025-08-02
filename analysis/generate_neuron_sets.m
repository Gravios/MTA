function neuronSets = generate_neuron_sets(trialPairNames, pairsCell, includeGuesses)
% generateNeuronSets  Build neuron identity sets that are common across trials
%
%   neuronSets = generateNeuronSets(trialPairNames, pairsCell)
%   neuronSets = generateNeuronSets(trialPairNames, pairsCell, includeGuesses)
%
% INPUTS
%   trialPairNames : M-by-M cell array; {i,j} = {'TrialA','TrialB'}
%   pairsCell      : M-by-M cell array; {i,j} = N-by-4 double
%                       [unitA, unitB, cellType, quality]
%   includeGuesses : logical scalar (default = false).  When false, rows
%                    whose quality == 0 are ignored.
%
% OUTPUT
%   neuronSets : 1-by-K cell array of tables.  Each table lists every
%                trial–unit instance that belongs to one neuron recorded
%                in ≥ 2 trials, together with its cell-type.
%
% EXAMPLE
%   load pairingData.mat   % defines trialPairNames & pairsCell
%   S = generateNeuronSets(trialPairNames, pairsCell);
% -------------------------------------------------------------------------

if nargin < 3,  includeGuesses = false;  end

%% 1. Initialise lookup containers
nodeKey = containers.Map('KeyType','char','ValueType','double'); % "Trial_Unit" → node index
trialList   = {};              % trial name  per node
unitIDList  = [];              % unit ID     per node
cellTypeList = [];             % cell-type   per node  (NaN until known)

src = [];   dst = [];          % edge list

[nRows, nCols] = size(trialPairNames);

%% 2. Build nodes and edges
for i = 1:nRows
    for j = 1:nCols
        pairNames = trialPairNames{i,j};
        if isempty(pairNames),  continue,  end
        
        T1 = pairNames{1};
        T2 = pairNames{2};
        P  = pairsCell{i,j};
        if isempty(P),          continue,  end
        
        for r = 1:size(P,1)
            if ~includeGuesses && P(r,4) == 0,  continue,  end  % skip low-quality
            
            u1 = P(r,1);  u2 = P(r,2);
            cType = P(r,3);   % assumed to be the same for both units
            
            key1 = sprintf('%s_%d',T1,u1);
            key2 = sprintf('%s_%d',T2,u2);
            
            % --- node 1
            if ~isKey(nodeKey,key1)
                idx = nodeKey.Count + 1;
                nodeKey(key1) = idx;
                trialList{idx,1}   = T1;
                unitIDList(idx,1)  = u1;
                cellTypeList(idx,1)= cType;
            else
                idx = nodeKey(key1);
                if ~isnan(cellTypeList(idx)) && cellTypeList(idx) ~= cType
                    warning('Cell-type conflict for %s (existing = %d, new = %d).', key1, cellTypeList(idx), cType);
                end
                cellTypeList(idx) = cType;   % overwrite NaN if necessary
            end
            
            % --- node 2
            if ~isKey(nodeKey,key2)
                idx = nodeKey.Count + 1;
                nodeKey(key2) = idx;
                trialList{idx,1}   = T2;
                unitIDList(idx,1)  = u2;
                cellTypeList(idx,1)= cType;
            else
                idx = nodeKey(key2);
                if ~isnan(cellTypeList(idx)) && cellTypeList(idx) ~= cType
                    warning('Cell-type conflict for %s (existing = %d, new = %d).', key2, cellTypeList(idx), cType);
                end
                cellTypeList(idx) = cType;
            end
            
            % --- add edge
            src(end+1,1) = nodeKey(key1);
            dst(end+1,1) = nodeKey(key2);
        end
    end
end

%% 3. Connected components
G  = graph(src,dst);
cc = conncomp(G);
K  = max(cc);

neuronSets = cell(1,K);
for k = 1:K
    members = find(cc == k);
    if numel(members) < 2,  continue,  end   % min set size = 2
    
    T = table(trialList(members), unitIDList(members), cellTypeList(members), ...
              'VariableNames', {'Trial','UnitID','CellType'});
    
    neuronSets{k} = T;
end
neuronSets = neuronSets(~cellfun('isempty',neuronSets));

fprintf('Generated %d neuron sets (size ≥ 2, guesses %s).\n', ...
        numel(neuronSets), ternary(includeGuesses,'included','excluded'));
end

% Small helper ------------------------------------------------------------
function out = ternary(cond,a,b)
out = b; if cond, out = a; end
end
