function Stc =  generate_stc_with_embedding_offset(Stc,trimLength)
%function Stc =  generate_stc_with_embedding_offset(Stc)
% 
% Trim edges of each state within an stc to prevent overlap in
% embedded analysis.
%
% Argin:
%
%     Stc: MTAStateCollection
%
%     trimLength: numeric, length in seconds to trim from ends of states
%

Stc.updateMode([Stc.mode,'_emb_trm',num2str(trimLength*1000)]);

for state = Stc.list_state_attrib('label'),
    Stc.states{Stc.gsi(state)} = [Stc{state{1}}]+[trimLength,-trimLength];
end

Stc.save(1);