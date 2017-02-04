function [Stc] = optimize_stc_transition_eduardo(Trial,varargin)

% Note that the values .3 for removing epochs and 1.2 for gaussian window
% were chossen by the lower errors rate (Zero gaps <6% and Behavior overlap <.5%)
% DEFARGS ------------------------------------------------------------------------------------------
[Stc,states,minChunkDuration] = DefaultArgs(varargin,{Trial.stc.copy,{'walk','rear','turn','pause','groom','sit'},[0.4,0.6,0.3,0.4,1.25,1.75]},true);

if ischar(Stc),
    Stc = Trial.load('stc',Stc);
end

%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
xyz = Trial.load('xyz');
nsts = numel(states);
Samples = round(minChunkDuration*xyz.sampleRate);
Binary=zeros(length(xyz.data),nsts);
BinaryConv=zeros(length(xyz.data),nsts);

for s = 1:nsts;
    State=[];
    State = Stc{states{s}}.data;

    % removing shorter behavior periods
    State(State(:,2)-State(:,1)<round(Samples(s)*.4),:)=[]; 
    % Making binary
    Tlength = length(xyz.data);        
    Binary= WithinRanges([1:Tlength],State);
    % Convolve with gaussian
    gauss=[];
    gauss = normpdf(linspace(0,2,round(Samples(s)*1.4)),1,2);
    gauss = gauss./sum(gauss);
    
    BinaryConv(:,s)=conv(double(Binary),gauss,'same');


end

[~,mind] = max(bsxfun(@rdivide, BinaryConv,max(BinaryConv)),[],2);

smat = zeros([size(xyz,1),nsts])';
smat([mind+(0:size(xyz,1)-1)'*nsts]) = 1;
smat = smat';
smat(sum(BinaryConv,2)<=0,:)=0;
smat = logical(smat);

%CorrectedBehavior = bsxfun(@rdivide, BinaryConv,max(BinaryConv))>.5;

for s=1:nsts;
    Stc.states{Stc.gsi(states{s})}.data = ThreshCross(smat(:,s),0.5,1);
end

Stc.updateMode([Stc.mode '_OPTE']);
Stc.save(1);

%---------------------------------------------------------------------------------------------------