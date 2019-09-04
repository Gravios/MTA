function diagnostic_lrzc_ssh(Trial)
Trial = MTATrial.validate(Trial);

disp('Starting: diagnostic_lrzc_ssh')
jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf'...
             ' -y ' Trial.path.project ' -l ' Trial.name ...
             ' diagnostic_lrzc_ssh_subp1']);
r1jid = [' -d afterok:' char(jid.readLine)]
