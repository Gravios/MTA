function diagnostic_lrzc_MatSubMitLRZ(Trial)
%function diagnostic_lrzc_MatSubMitLRZ(Trial)
%testing matlab function submission to the lrz clusters

Trial = MTATrial.validate(Trial);
cd(Trial.path.data);


%1. 
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp1')
jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf -l ' Trial.name ...
             ' diagnostic_lrzc_MatSubMitLRZ_subp1']);
r1jid = [' -d afterok:' char(jid.readLine)]

%2. 
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp2')
system(['MatSubmitLRZ --config lrzc_hugemem.conf '...
         r1jid ' -l ' Trial.name ...
        ' diagnostic_lrzc_MatSubMitLRZ_subp2'])
%3.
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp3')
jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf '...
                   r1jid ' -l ' Trial.name ...
                  ' diagnostic_lrzc_MatSubMitLRZ_subp3']);
r3jid = [' -d afterok:' char(jid.readLine)]

%4.
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp4')
system(['MatSubmitLRZ --config lrzc_hugemem.conf' ...
         r3jid ' -l ' Trial.name ...
        ' diagnostic_lrzc_MatSubMitLRZ_subp4'])
