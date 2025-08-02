function [Stc] = optimize_stc_transition_anton(Trial,varargin)
% This function fixes the output of the Neural network data removing
% shorting period and filling gaps.
% Trial should be MTA object with loaded stc data
% minChunkDuration: Threshold base in distribution hand labeled behavior 
% the order is Walk,Rear,Trun,Pause,Groom,Sit.

% DEFARGS ------------------------------------------------------------------------------------------
[Stc,states,minChunkDuration] = DefaultArgs(varargin,{Trial.stc.copy,{'walk','rear','turn','pause','groom','sit'},[0.4,0.6,0.3,0.4,1.25,1.75]},true);

if ischar(Stc),
    Stc = Trial.load('stc',Stc);
end

%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
Chunks = [];
for s = 1:6,
    state = Stc{states{s}};
    Chunks = [Chunks; [state.data diff(state.data,1,2), s*ones([size(state,1),1])]];
end

Chunks = sortrows(Chunks,1);
ChunkNeighborsSame = [0; Chunks(3:end,4) ==Chunks(1:end-2,4);0];
ChunksAreShort = zeros(size(Chunks,1),1);

for s = 1:6,
    state = Stc{states{s}};    
    ChunkID = Chunks(:,4) == s;
    ChunksAreShort(ChunkID) = (Chunks(ChunkID,3)<minChunkDuration(s)*state.sampleRate);
end

ChunksFixed = Chunks;
ChunksFixed(ChunksAreShort & ~ChunkNeighborsSame,4) = 0;
ChunksFixed(ChunksAreShort &  ChunkNeighborsSame,4) = ChunksFixed( find(ChunksAreShort & ChunkNeighborsSame)+1,4);

for s=1:6;
    Stc.states{Stc.gsi(states{s})}.data = ChunksFixed(ChunksFixed(:,4)==s,1:2);
end

Stc.updateMode([Stc.mode '_OPTA']);
Stc.save(1);
%---------------------------------------------------------------------------------------------------