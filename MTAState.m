classdef MTAState
%MTAState(key,lable,state) - container for behavioral state information
%
%key - string: keyboard character associated with state
%
%label - string: name describing state (e.g. 'rear' or 'walk')
%
%state - numericArray<-double: (event,[time_onset,time_offset]) @xyzSampleRate
%

    properties (SetAccess = public)

        %key - string: keyboard character associated with state
        key

        %label - string: name describing state (e.g. 'rear' or 'walk')
        label

        %state - numericArray<-double: (event,[time_onset,time_offset]) @xyzSampleRate
        state

    end

    methods

        function State = MTAState(key,label,state)
            State.key = key;
            State.label = label;
            State.state = state;
        end

    end

end