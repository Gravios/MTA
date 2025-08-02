function Par = load(Par,Session)
% function Par = load(Par)
% load(Par)
% loads the specified par file and returns a structure with these elements:
%
% .FileName      -> name of file loaded from
% .nChannels     -> number of total channels
% .nBits         -> number of bits of the file
% .SampleTime    -> time, in microseconds, of 1 sample (ie 1e6 / sample rate)
% .HiPassFreq    -> High pass filter frequency
% .nElecGps      -> number of electrodes (i.e. electrode groups)
% .ElecGp        -> a cell array giving the channels in the electrodes
%                    e.g. if .ElectrodeGroup{3} = [2 3 4 5], electrode 3
%                    is a tetrode for channels 2 3 4 and 5. 
% channel numbers here are from 0. be carefull.
%

Par.data = load_xml(Par.fpath());
Par.parent = Session;
