function Data = create(Data,Session,varargin)
%create(Data,Session,varargin)
%[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
Data.load(Session,varargin{:});
end
