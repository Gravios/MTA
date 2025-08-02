function diagnostic_lrzc_MatSubMitLRZ_parfor_subp1(Trial)
Trial = MTATrial.validate(Trial);
disp(['hello from ',mfilename,Trial.filebase,' was my input'])
disp(['Opening parpool with: ' num2str(nWorkers) ' workers'])
pp = parpool(4);
parfor f = 1:5,
    pause(300);
    disp(['iter: ',num2str(f),' in ',mfilename,'... complete'])
end
disp(['Closing parpool'])
delete(pp)
disp(['exiting: ',mfilename])
