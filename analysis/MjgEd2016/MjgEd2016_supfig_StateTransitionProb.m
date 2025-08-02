




transWindow = repmat({0.2},1,numSessions);
tsts = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(tsts);
stpa = {};
for t = 1:nsts,
    for o = 1:nsts,
        if t~=o,
            %   SELECT state to state transitions
            stpa{t,o} = cell2mat(cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{tsts{t},tsts{o}},tw,x),Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz)');
        end
    end
end

stf = cell2mat(cf(@length,stpa));
round(bsxfun(@rdivide,stf,sum(stf)),2)
round(bsxfun(@rdivide,stf,sum(stf,2)),2)
stp = round(stf./sum(stf(:)),4);

%         0    0.0224    0.0391   0.2591  0.0022  0.0088
%    0.0099         0    0.0177   0.0188  0.0011       0
%    0.0893    0.0011         0   0.0757  0.0011  0.0004
%    0.2157     0.019    0.1012        0  0.0237  0.0291
%    0.0041    0.0015    0.0019   0.0201       0       0
%     0.011    0.0004    0.0041   0.0207  0.0006       0



for t = 1:nsts
    for o = 1:nsts
        
    end
end
