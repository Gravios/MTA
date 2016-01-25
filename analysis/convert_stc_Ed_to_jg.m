function convert_stc_Ed_to_jg(Trial)

if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end

Trial.load('stc','manual1');

Trial.stc.states{Trial.stc.gsi('pause')}.key = 'p';
Trial.stc.states{ Trial.stc.gsi('walk')}.key = 'w';
Trial.stc.states{ Trial.stc.gsi('rear')}.key = 'r';
Trial.stc.states{ Trial.stc.gsi('turn')}.key = 'n';
Trial.stc.states{Trial.stc.gsi('shake')}.key = 'k';
Trial.stc.states{  Trial.stc.gsi('sit')}.key = 's';
Trial.stc.states{Trial.stc.gsi('groom')}.key = 'm';

Trial.stc.updateMode('hand_labeled_rev1_Ed');
Trial.stc.save(1);