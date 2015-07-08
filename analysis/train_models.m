%
% For now the only well labeled session is 'jg05-20120317'
% Most up to date Stc is 'hand_labeled_rev1'
%

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev1');

bhv_lgr(Trial,true);
bhv_lda(Trial,true);

[Stc,d_state] = bhv_lgr(Trial);
Stc.save(1);
[Stc,d_state] = bhv_lda(Trial);
Stc.save(1);