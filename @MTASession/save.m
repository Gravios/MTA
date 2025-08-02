function save(Session)
%save(Session)
save(fullfile(Session.spath, [Session.filebase '.ses.mat']),'Session','-v7.3');
end

