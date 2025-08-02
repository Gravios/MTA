function diagnostic_lrzc_MatSubMitLRZ_args(Trial)
%function diagnostic_lrzc_MatSubMitLRZ(Trial)
%testing matlab function submission to the lrz clusters

Trial = MTATrial.validate(Trial);

%1. 
disp('Starting: diagnostic_lrzc_MatSubMitLRZ_args_subp1')
jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf'...
             ' -y ' Trial.path.project ' -l ' Trial.name ...
             ' diagnostic_lrzc_MatSubMitLRZ_subp1 69 \''bob\''']);

%system('MatSubmitLRZ -v --config lrzc_hugemem.conf -y /home/hpc/pr84qa/di68tor/data/project/general -l jg05-20120317 diagnostic_lrzc_MatSubMitLRZ_args_subp1 69 \''bob\''')

%r1jid = [' -d afterok:' char(jid.readLine)]

% $$$ %2. 
% $$$ disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp2')
% $$$ system(['MatSubmitLRZ --config lrzc_hugemem.conf '...
% $$$         ' -y ' Trial.path.project r1jid ' -l ' Trial.name ...
% $$$         ' diagnostic_lrzc_MatSubMitLRZ_subp2'])
% $$$ %3.
% $$$ disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp3')
% $$$ jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf '...
% $$$                   ' -y ' Trial.path.project r1jid ' -l ' Trial.name ...
% $$$                   ' diagnostic_lrzc_MatSubMitLRZ_subp3']);
% $$$ r3jid = [' -d afterok:' char(jid.readLine)]
% $$$ 
% $$$ %4.
% $$$ disp('Starting: diagnostic_lrzc_MatSubMitLRZ_subp4')
% $$$ system(['MatSubmitLRZ --config lrzc_hugemem.conf' ...
% $$$         ' -y ' Trial.path.project r3jid ' -l ' Trial.name ...
% $$$         ' diagnostic_lrzc_MatSubMitLRZ_subp4'])
