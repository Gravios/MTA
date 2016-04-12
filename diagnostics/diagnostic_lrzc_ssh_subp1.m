function diagnostic_lrzc_ssh_subp1(Trial)

Trial = MTATrial.validate(Trial);
disp('Starting: diagnostic_lrzc_ssh_subp1')
disp(['hello from ',mfilename, ' ',Trial.filebase,' was my input'])
popen(['ssh di68tor@lxlogin1.lrz.de \"MatSubmitLRZ --config lrzc_hugemem.conf'...
             ' -y ' Trial.path.data ' -l ' Trial.name ...
             ' diagnostic_lrzc_ssh_subp1\"']);
r1jid = [' -d afterok:' char(jid.readLine)]
system(['ssh di68tor@lxlogin1.lrz.de \"sbatch /naslx/projects/pr84qa/di68tor/code/sbatch/diagnostic_lrzc_ssh_subp1_jg05-20120317/diagnostic_lrzc_ssh_subp1_jg05-20120317.cmd\"'])



