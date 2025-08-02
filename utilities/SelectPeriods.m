%function [y, ind] = SelectPeriods(x,Periods,SigType,WhereFlag, ifSquash)
% selects from the signal (time series  - 'd', continuous - 'c', 'p' - periods)
% values in/out the Periods ranges (should be the same sampl.rate
% WhereFlag defines in(1) or out(0) (default - 0 = in)
% ifSquash is 1 if you apply it for discrete signal and want to squash the
% gaps to match the sample indexes of continuous signal
%   output: y - new signal, ind - indexes of the old signal
function [y, ind] = SelectPeriods(x,Periods,SigType,varargin)
[WhereFlag, ifSquash] = DefaultArgs(varargin,{1,0});

% IF periods is a string assume it is a filename and load the file
if isstr(Periods)
    Periods = load(Periods);
end

% IF no periods are specified return x
if isempty(Periods)
    y = x;
    ind = [1:length(x)];
    return
end

% REMOVE any periods contaning zero(s)
Periods(find(Periods(:)==0))=1;

nPeriods  = size(Periods,1);

if ~iscell(x)
    nChannels = min(size(x));
    if size(x,1)==nChannels
        x=x';
    end
    nTimeBins = max(size(x));
    
    if (nargin<3 | isempty(SigType))
        if (nTimeBins < max(Periods(:)))
            SigType = 'd';
        else
            SigType = 'c';
            if (nargout >2)
                error('too many output parameters');
                exit;
            end
        end
    end
    if (SigType =='d')
        Channels = size(x,1);
    end
else
    celli=size(x,1);
    cellj=size(x,2);
end


if (nargin <4 | isempty(WhereFlag) )
    WhereFlag =1;
end

switch SigType
    case 'd'
        
        if ~iscell(x)
            [y ind] = SelPerDiscr(x,Periods,WhereFlag,ifSquash);
        else
            y = cell(celli,cellj);
            ind = cell(celli,cellj);
            for ii=1:celli
                for jj=1:cellj
                    [y{ii,jj} ind{ii,jj}]= SelPerDiscr(x{ii,jj},Periods,WhereFlag,ifSquash);
                end
            end
        end
        
        
        
        
    case 'c'
        y=[];
        ind=[];
        
        if WhereFlag
            for p=1:nPeriods
                y =   cat(1,  y, x(Periods(p,1):Periods(p,2),:,:,:,:));
                ind = cat(1,ind, [Periods(p,1):Periods(p,2)]');
            end
        else
            if Periods(1,1)>2
                y =   cat(1,  y, x(1:Periods(1,1)-1,:,:,:,:));
                ind = cat(1,ind, [1:Periods(1,1)-1]');
            end
            for p=1:nPeriods-1
                y =   cat(1,  y, x(Periods(p,2)+1:Periods(p+1,1)-1,:,:,:,:));
                ind = cat(1,ind, [Periods(p,2)+1:Periods(p+1,1)-1]');
            end
            
            y =   cat(1,  y, x(Periods(nPeriods,2)+1:end,:,:,:,:));
            ind = cat(1,ind, [Periods(nPeriods,2)+1:size(x,1)]');
            
        end
        
    case 'p' %this case is for periods
        if size(x,2)~=2
            error('x has to be two column matrix - periods of [beg end]');
        end
        if WhereFlag 
            x_Periods = IntersectRanges(x,Periods);
            % x_Periods should be inclusive Periods
            
            [yl indl] = SelPerDiscr(x_Periods(:,1),Periods,1,ifSquash);
            [yr indr] = SelPerDiscr(x_Periods(:,2),Periods,1,ifSquash);
            if ~isempty(setdiff(indl,indr))
                error('smth went wrong');
            end
            y = [yl yr];
            
        else
            error('this function is not implemented');
        end
       
    otherwise
        error('Parameter "SigType" must be either "c" or "d" or "p"!')        
end

