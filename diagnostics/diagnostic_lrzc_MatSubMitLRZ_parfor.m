function diagnostic_lrzc_MatSubMitLRZ_parfor(Trial)
%function diagnostic_lrzc_MatSubMitLRZ(Trial)
%testing matlab function submission to the lrz clusters

Trial = MTATrial.validate(Trial);

%1. SERIAL
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_parfor_subp1')
jid = popen(['MatSubmitLRZ --config lrzc_serialS.conf -l ' Trial.name ...
             ' diagnostic_lrzc_MatSubMitLRZ_subp1']);
r1jid = [' -d afterok:' char(jid.readLine)]

%2. MPP1
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_parfor_subp1')
jid = popen(['MatSubmitLRZ --config lrzc_mpp1S.conf -l ' Trial.name ...
             ' diagnostic_lrzc_MatSubMitLRZ_subp1']);
r1jid = [' -d afterok:' char(jid.readLine)]

%2. MYR1
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_parfor_subp1')
jid = popen(['MatSubmitLRZ --config lrzc_myr1.conf -l ' Trial.name ...
             ' diagnostic_lrzc_MatSubMitLRZ_subp1']);
r1jid = [' -d afterok:' char(jid.readLine)]
