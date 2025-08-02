function [StcCor,tempDState,tempPState] = reassign_state_by_duration(StcCor,key,defaultState,durationThreshold,tempDState,tempPState,logicFun)
pd = [];
states = StcCor.list_state_attrib;
pthresh = log10(durationThreshold.*StcCor{key}.sampleRate);
for rp = StcCor{key}.data',

    % collect period durations
    pd(end+1) = log10(abs(diff(rp)));

    if logicFun(pd(end),pthresh),
        pind = rp(1):rp(2);

        %???
        tempDState(pind,2) = 0;        
        % promote next best state 
        tempPState(pind,:) = circshift(tempPState(pind,:),-1,2);

        if isempty(defaultState), 
            % select next best state             
            msts = mode(tempPState(pind,1));
            if msts == 2; 
                tps(pind,:) = circshift(tempPState(pind,:),-1,2);
                msts = mode(tempPState(pind,1));
            end                    
        else % relabel state with a provided default state            
            msts = StcCor.gsi(defaultState);
        end
        
        
        % reassign state
        StcCor.states{StcCor.gsi(states{msts})}.data = ...
            [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
    end
end

% Resort and clean overlaps 
for sts = StcCor.states,
    sts{1}.clean; 
end

% Delete the reassigned periods within original state
StcCor.states{StcCor.gsi(key)}.data(logicFun(pd,pthresh),:) = [];
