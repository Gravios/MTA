function diagnostic_lrzc_MatSubMitLRZ_args_subp1(Trial,num,str)
Trial = MTATrial.validate(Trial);
disp(['hello from ',mfilename,strjoin({Trial.filebase,num,str},', '),' was my input'])
disp(['exiting: ',mfilename])