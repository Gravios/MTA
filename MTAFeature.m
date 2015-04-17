classdef MTAFeature
%MTAFeature(lable,expression,ifsave) - container for behavioral state information
%
%key - string: keyboard character associated with state
%
%label - string: name describing state (e.g. 'rear' or 'walk')
%
%expression - string: matlab expression defining the feature;
%
%ifsave - logical: save feature vector if true
%
%feature - numericArray: feature timeseries 

    properties (SetAccess = public)

        %key - string: keyboard character associated with state
        key

        %label - string: name describing state (e.g. 'rear' or 'walk')
        label
        
        plot_type

        %expression - string: matlab expression defining the feature;
        expression 
        
        %ifsave - logical: save feature vector if true
        ifsave
        
        %feature - numericArray: feature timeseries
        feature        

    end

    methods

        function Feature = MTAFeature(label,expression,ifsave,plot_type)
            Feature.label = label;
            Feature.expression = expression;
            Feature.ifsave     = ifsave;
            Feature.feature    = [];
            Feature.plot_type  = plot_type;
        end

        function Feature = compute(Feature,Session)            
            Feature.feature = eval(Feature.expression);
            if isa(Feature.feature,'MTAData'),
                Feature.feature = Feature.feature.data;
            end
            
        end

    end

end