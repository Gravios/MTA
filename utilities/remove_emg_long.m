function varargout = remove_emg_long(Session, fileExt, data, varargin)
% function [data, Ws, As, EMG_au, AW, armodel] = EMG_rm_long(x,[
%                                       denoise_frequency_lowerbound,
%                                       silence_periods,included_periods,
%                                       SamplingRate, 
%                                       high_pass_freq, EMG_thrd, if_rm_mean, 
%                                       armodel, cmp_method, down_sample,
%                                       isave, save_range, FileName, savedir,
%                                       save_together,use_wb,EMG_shank,EMGFileName])
% function to remove the EMG noise in a long period. 
% EMG detected by high frequency (>high_pass_freq) correlated activity
% among channels. 
% Inputs: 
%   x: data, nt x nch.
%   Optional: 
%       denoise_frequency_lowerbound: keep slower signal until this frequency
%       silence_periods: whether one would like to remove noise only for
%           non_silenced periods
%       included_periods: the periods to be included.
%       SamplingRate: in Hz, default: 1000 Hz.
%       removeLinenoiseFlag: if remove line noise. default: true
%       line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.
%       high_pass_freq: beginning of frequency to detect the muscle tone
%                       default: 100 Hz.
%       EMG_thrd: the thresholded (binary) data of the high EMG periods. 
%               Precomputed by the function of EMG_Cluster.m.
%       if_rm_mean: if return the decenterized signal, default: true.
%       armodel: the armodel or the name of the armodel file. if not given,
%               we woud use whiten signal to compute a corresponding
%               ARmodel. 
%       cmp_method: methods to compute the EMG components. Whiten('w') or
%                   Highpass&Whiten('hw', a bit more stable). defualt: 'hw'
%       down_sample:down sample the data to compute the periods.
%                   default:3.
%       numOfIC: adjusted according to number of channels.
%       isave: if save to a file. defualt: false (when the entry is empty)
%              if the entry is not empty, the file is saved to a file
%              given by "Filename".
%       save_range: the Period in the original whole dataset, used to write
%                   the cleaned data into a larger binary file. 
%                   formate: cell(2,1): 
%                           save_range{1}: [Start_Sample  End_Sample;
%                                           Start_Channel End_Channel];
%                           save_range{2}: [total_channels, total_samples],
%                   if not given, save as defualt: [1 nt;1 nch]. 
%       FileName: defualt:  current_dictionary_name.EMG_rm.mat and
%                           current_dictionary_name.lfpc.
%       savedir:    the directory that one is gonna save the data. the
%                   default is the current directory.
%       save_together: if save ICs from all the chunks together. 
%                       default: true
%       use_wb: whether to use the wide band signal when the algorithm
%               can't find a proper flat component at high frequency band
%               defualt: true NB: please check the report fig to see if you
%               really need this! 
%       EMG_shank: the shank currently used. 0: no shank number is given
%       HP:
%       Shk01:
%       EMGfile
%       
% Outputs:
%   x: EMG removed signal.
%   Ws: unmixing vector of the EMG component.
%   As: loading of the EMG component.
%   EMG_au: normalized EMG component in arbitrary unit.
%   AW: everything about the ICA 
%   armodel: the armodel used. (use the same armodel for all the chunks)
%
%   Note:   function without output arguments saves the EMG component and
%           the binary lfp data to files. For large data I use memmapfile.
% 
% Related functions: 
% EMG_Cluster.m, SREaffineV.m, EMG_rm_main.m.
% This function is a part of the EMG_removing toolbox.
%  
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 01.12.2019.

% DEFARGS ---------------------------------------------------------------------------------------
defargs = struct(                                                                             ...
    'denoiseFrequencyLowerbound', 0,                                                          ...
               'silencePeriods', false,                                                       ...
              'included_periods', [],                                                         ...
                'dataSampleRate', 1250,                                                       ...
                  'removeLinenoiseFlag', true,                                                ...
                     'line_thrd', 2,                                                          ...
      'denoiseFrequencyHighpass', 100,                                                        ...
                      'EMG_thrd', [],                                                         ...
                    'if_rm_mean', true,                                                       ...
                       'armodel', [],                                                         ...
                    'cmp_method','hw',                                                        ...
                   'down_sample', 3,                                                          ...
                       'numOfIC', 0,                                                          ...
                         'isave', false,                                                      ...
                    'save_range', [],                                                         ...
                 'save_together', true,                                                       ...
                        'use_wb', false,                                                      ...
                     'probeId', 0,                                                            ...
                            'HP', [],                                                         ...
                         'Shk01', [],                                                         ...
                   'filePathEMG', [],                                                         ...
                'filePathOutput', []                                                          ...
);
[denoiseFrequencyLowerbound,silencePeriods, included_periods, dataSampleRate, removeLinenoiseFlag,  ...
 line_thrd,denoiseFrequencyHighpass, EMG_thrd, if_rm_mean,armodel,cmp_method,down_sample,     ...
 numOfIC,isave,save_range,save_together,use_wb,probeId,HP,Shk01,filePathEMG,filePathOutput] = ... 
    DefaultArgs(varargin,defargs,'--struct');
%------------------------------------------------------------------------------------------------


[nt, nch] = size(data);
if isempty(included_periods)
    included_periods = 1:nt;
else
    sds = StartEnd1d(included_periods);
    included_periods = find(included_periods);
end

if nt < nch 
    istr = true;
    data = data';
    [nt, nch] = size(data);
else % why? check later
    istr = false;
end

assert( ~isempty(EMG_thrd), 'MTA:utilities:remove_emg_long:HighEmgDetectionFileNotFound.');

if ~numOfIC
    numOfIC = max(fix(.7*nch), min(10,nch));
end

selectedprd = EMG_thrd;

%% INITIALIZE the OUTPUTS
Ws= zeros(1,nch);
As= zeros(nch,1);
EMG_au= zeros(nt,1);
%% PREPARE DATA

mx = mean(data);
data = bsxfun(@minus, data, mx);

% Whitening data
if isempty(armodel)
    armodel = dir(fullfile(savedir,[Session.name,'.whitensignal.ARmodel.',fileExt,'.mat']));
    if isempty(armodel)
        [wx,armodel]=whitensignal(data,[],[],[],1);
        tmp_ar = armodel;
    else
        armodel = load(armodel(1).name);
        wx=whitensignal(data,[],[],armodel.ARmodel);
        tmp_ar = armodel.ARmode;
    end
elseif isstr(armodel)
    armodel = load(armodel);
    tmp_ar = armodel.ARmodel;
    wx=whitensignal(data,[],[],tmp_ar);
elseif isfield(armodel,'ARmodel')
    tmp_ar = armodel.ARmodel;
    wx=whitensignal(data,[],[],tmp_ar);
else 
    try
        tmp_ar = armodel;
        wx=whitensignal(data,[],[],armodel);
    catch 
        warning('Something wrong with the AR model you give, try whiten without given model.')
        [wx, armodel]=whitensignal(data,[],[],[],1);
        tmp_ar = armodel;
    end
end


AW.armodel = tmp_ar;
% opf_A = @(x)(bsxfun(@rdivide,x,std(x)));
if isempty(HP)
    chmap = 1:nch;
    if isempty(save_range)
        HP=chmap;
    else
        HP = save_range{1}(2,1):save_range{1}(2,2);
    end
else
    chmap=HP;
end
chmap = [1:4,6:64];

opf_a  = @(x) bsxfun(@rdivide,x,sqrt(sum(x.^2)));
opf_Aa = @(x) bsxfun(@rdivide,x,remove_emg_SREaffineV(chmap,x,Shk01));
opf_A  = @(x) opf_Aa(opf_a(x));
% SREaffineV use affine to fit, and compute the variance accordinglty.
% Accounting for the linear leaking from other areas.  
%% REMOVE LINE-NOISE
% note: line noise is still removed for any periods. 
if removeLinenoiseFlag
    [A_rm_line,W_rm_line,A_line,W_line,power_ratio,thrd] = remove_emg_linenoise(wx,line_thrd,dataSampleRate);
    AW.A_rm_line = A_rm_line;
    AW.W_rm_line = W_rm_line;
    AW.A_line = A_line;
    AW.W_line = W_line;
    AW.power_ratio = power_ratio;
    AW.thrd = thrd;
    if ~isempty(A_rm_line)
        signals = data * W_rm_line';
        for n = 1:size(signals,2)
            tmp_id = find(sign(signals(1:(end-1),n)).*sign(signals(2:end,n))<=0, 1,'first');
            signals(1:tmp_id,n) = 0;
            tmp_id = find(sign(signals(1:(end-1),n)).*sign(signals(2:end,n))<=0, 1,'last');
            signals(tmp_id:end,n) = 0;
        end
        data = data - signals*A_rm_line';
        wx = wx - wx*W_rm_line'*A_rm_line'; % won't affect the results. 
        fprintf('\n Line noise component removed...\n')
    end
end

%% EMG COMPONENTS AND ACTIVITIES
flatness_threshold = 2*nch;

% What is the reasoning behind this threshold
if sum(selectedprd)>=(10*dataSampleRate)
    remove_cmp = true;
    AW.usewb = false;
    % Components from the high frequency.
    switch lower(cmp_method) %
      case 'hw'
        % defualt
        wx = wx(:,chmap);
        hx = ButFilter(data(:,chmap,4,[denoiseFrequencyHighpass,1000]/(dataSampleRate/2),'bandpass');
        [~,sind] = sort(sum(hx(selectedprd,:).^2));
        [Ah, Wh] = fastica(hx(selectedprd,sind)', 'numOfIC', numOfIC,'verbose','off');

        [~,rod] = sort(abs(sum(opf_A(Ah))),'descend');
        % Components from all
        %ICA computed on the highpass filtered data using the selected periods
        
        [Ax, Wx] = fastica(Wh(rod,:)*wx(included_periods(1:down_sample:end),sind)','verbose','off');
        
        A = Ah(:,rod)*Ax;
        W = Wx*Wh(rod,:);
        flatness = abs(sum(opf_A(A)));
        if (sum(abs(sum(opf_A(A)))>nch)>1) || sum(abs(sum(opf_A(A)))>flatness_threshold)<1 
            remove_cmp = false;
            fprintf('recompute...\n')
            wx=WhitenSignal(wx,[],[],tmp_ar);
            nn = 1;
            
            while  ((sum(flatness>nch)>1) || sum(flatness>flatness_threshold)<1)&&(nn<5)
                % too many flat components or none of the components are flat enough.
                if ((sum(abs(sum(opf_A(Ah)))>nch)>1) || sum(abs(sum(opf_A(Ah)))>flatness_threshold)<1)
                    [Ah, Wh] = fastica(hx(selectedprd,:)','verbose','off'); % 'numOfIC', numOfIC,
                end
                [Ax, Wx] = fastica(Wh*wx(included_periods(nn:down_sample:end),:)','verbose','off');
                A = Ah*Ax;
                W = Wx*Wh;
                nn = nn+1;
                fprintf('\r%d in %d...\n',nn,5)
                flatness = abs(sum(opf_A(A)));
            end
            if ((sum(abs(sum(opf_A(Ah)))>nch)<2) && sum(flatness>flatness_threshold)>0)
                % too many components usually indicates there's no typical EMG noise. 
                remove_cmp = true;
            end
            
            fprintf('\n\n')
        end
        if (sum(flatness>nch)<1) && use_wb
            [Awb, Wwb] = fastica( data(included_periods(1:down_sample:end),:)', 'verbose','off');% se
            AW.Awb = Awb;
            AW.Wwb = Wwb;
            if sum(abs(sum(opf_A(Awb)))>(nch*.9))>0
                A = Awb;
                W = Wwb;
                AW.usewb = true;
                fprintf('\r Using the non-whitened data.\n')
                remove_cmp = true;
            end
        end
        AW.Ah = Ah;
        AW.Ax = Ax;
        AW.Wh = Wh;
        AW.Wx = Wx;
      case 'w' % Use Whiten alone, pretty loose
        [A, W] = fastica(wx(included_periods(1:down_sample:end),:)', 'numOfIC', numOfIC,'verbose','off');% selectedprd
        flatness = abs(sum(opf_A(A)));
        if (sum(flatness>nch)<1)
            remove_cmp = false;
            if use_wb
                [Awb, Wwb] = fastica(data(included_periods(1:down_sample:end),:)', 'verbose','off');% se
                AW.Awb = Awb;
                AW.Wwb = Wwb;
                if sum(abs(sum(opf_A(Awb)))>(nch*.9))>0
                    A = Awb;
                    W = Wwb;
                    AW.usewb = true;
                    fprintf('\r Using the non-whitened data.\n')
                    remove_cmp = true;
                end
            end
        end
      otherwise
        fprintf('Please using hw or w. ')
    end
    
    flatness = abs(sum(opf_A(A)));
    [~, EMG_comp] = max(flatness);
    As = A(:,EMG_comp);
    Ws = W(EMG_comp,:);
    if denoiseFrequencyLowerbound>0
        EMG_au = ButFilter((data * Ws'), 4, denoiseFrequencyLowerbound/(dataSampleRate/2),'high');%.*selectedprd;
    else
        EMG_au = data * Ws';%.*selectedprd;
    end
    % (:,n) is to see the behavior of the highpassed resulted EMG.
    % not need in the final varsion.
    
    AW.A = A;
    AW.W = W;
    
    if remove_cmp
        %% SMOOTHING THE BOUNDRIES OF SILENCED PERIODS.
        cross_searching_ranges = 50;
        if silencePeriods
            last_end = 1;
            for k = 1:size(sds,1)
                % for the starting points
                if sds(k,1)<cross_searching_ranges
                    [~,id] = remove_emg_find_crossing(abs(EMG_au(1:cross_searching_ranges)),1);
                    EMG_au(1:id) = 0;
                else
                    tmp = (sds(k,1)-cross_searching_ranges):sds(k,1);
                    [~,id] = remove_emg_find_crossing(abs(EMG_au(tmp)),-1);
                    EMG_au(last_end:tmp(id)) = 0;
                end
                % for the ends
                if sds(k,2)<(nt-cross_searching_ranges)
                    tmp = [0:cross_searching_ranges]+sds(k,2);
                    [~,id] = remove_emg_find_crossing(abs(EMG_au(tmp)),1);
                    last_end = tmp(id)+1;
                end
            end
            if last_end>1
                EMG_au(last_end:end) = 0;
            end
        end
        %% SMOOTHING THE ENDS OF EMG TRACES.
        
        % first_cross
        if abs(EMG_au(1))<(median(abs(EMG_au))*1e-3)
            first_cross = 1;
        else
            first_cross = find(sign(EMG_au(1:(end-1)).*EMG_au(2:end))<0,1,'first');
        end
        if first_cross>1250
            first_cross = find(EMG_au<=(median(abs(EMG_au))*1e-2),1,'first');
        end
        EMG_au(1:first_cross)=0;
        AW.zero_first = first_cross;% zeros line at first.
        
        % last_cross
        if abs(EMG_au(end))<(median(abs(EMG_au))*1e-3)
            last_cross = nt;
        else
            last_cross = find(sign(EMG_au(1:(end-1)).*EMG_au(2:end))<0,1,'last');
        end
        if (nt-last_cross)>1250
            last_cross = find(EMG_au<=(median(abs(EMG_au))*1e-2),1,'last');
        end
        EMG_au((last_cross+1):end)=0;
        AW.zero_last = nt - last_cross+1;% zeros line at last.
        
        data = data - EMG_au*As';% (:,end)
        
    else
        EMG_au = EMG_au*0;
        fprintf('EMG_component is not removed')
    end
    if ~if_rm_mean
        data = bsxfun(@plus, data, mx);
    end
    if istr % transpose back.
        data = data';
    end
end
%% OUTPUT OR SAVE DATA

par.denoiseFrequencyLowerbound = denoiseFrequencyLowerbound;
par.silencePeriods =silencePeriods;
par.included_periods =included_periods;
par.dataSampleRate =dataSampleRate;
par.removeLinenoiseFlag =removeLinenoiseFlag;
par.line_thrd =line_thrd;
par.denoiseFrequencyHighpass =denoiseFrequencyHighpass;
par.EMG_thrd =EMG_thrd;
par.if_rm_mean=if_rm_mean;
par.armodel=armodel;
par.cmp_method=cmp_method;
par.down_sample=down_sample;
par.numOfIC=numOfIC;
par.isave=isave;
par.save_range=save_range;
par.save_together=save_together;
par.use_wb=use_wb;
par.probeId=probeId;
par.HP=HP;
par.Shk01=Shk01;
par.filePathEMG=filePathEMG;
par.remove_cmp = remove_cmp;
par.flatness = flatness;
scaling_factor = sign(sum(As))*sqrt(As'*As);
if nargout>1
    varargout{1} = data;
    varargout{2} = Ws;
    varargout{3} = As;
    varargout{4} = EMG_au;
    varargout{5} = AW;
    varargout{6} = tmp_ar;
    varargout{7} = scaling_factor;
    varargout{8} = flatness;
    varargout{9} = par;
else 
    isave = true;
end

% IGNORE for now - isave will remain false
% function without output arguments saves the EMG component and the binary
% lfp data to files. 
if isempty(save_range)
    save_range =    [1 nt;...
                    1 nch];
end

% SAVE intermittently
if isave
% $$$     if isempty(FileName)
% $$$ % $$$         if ~isempty(savedir)
% $$$ % $$$             cwd = savedir;
% $$$ % $$$         else
% $$$ % $$$             cwd = pwd;
% $$$ % $$$         end
% $$$ % $$$         FileName =  cwd((find(cwd=='/',1,'last')+1):end);
% $$$         
% $$$         filePathOutput = fullfile(savedir,[Session.name,'.',fileExt,'d']);
% $$$         
% $$$         if isempty(filePathEMG)
% $$$             filePathEMG = sprintf('%s%s.sh%d.emg',savedir,FileName,probeId);
% $$$         end
% $$$         
% $$$         if exist(filePathOutput,'file')
% $$$             warning(sprintf('%s already exist! Please check!\n Now saving to the %s.new files.', FileName, FileName))
% $$$             FileName = [FileName,'.new'];
% $$$             filePathOutput = sprintf('%s%s.lfpd',savedir,FileName);
% $$$         end
% $$$         if ~save_together
% $$$             G_par.HPs = HP;
% $$$             G_par.Shk01=Shk01;
% $$$             denoise_shank=probeId;
% $$$             save(sprintf('%s.EMG_rm.t%d-%d.ch%d-%d.mat',FileName,save_range(1,:),save_range(2,:)), 'AW','EMG_au','armodel','Ws','As','armodel','selectedprd','denoise_shank','filePathOutput','scaling_factor','G_par','par','flatness')
% $$$         end
% $$$         fileID = fopen(filePathOutput,'w');
% $$$         fwrite(fileID, int16(data'),'int16');
% $$$         fclose(fileID);
% $$$         
% $$$         fileID = fopen(filePathEMG,'w');
% $$$         tmp_EMG_au = scaling_factor*EMG_au;
% $$$         fwrite(fileID, int16(tmp_EMG_au'),'int16');
% $$$         fclose(fileID);
% $$$     else
% SAVE parameters
        if ~save_together
            G_par.HPs     = HP;
            G_par.Shk01   = Shk01;
            save(fullfile(savedir,                                                                                    ...
                          [Session.name,sprintf('.EMG_rm.t%d-%d.ch%d-%d.mat',save_range{1}(1,:),save_range{1}(2,:))]),...
                 'AW','EMG_au','armodel','Ws','As','armodel','selectedprd','probeId',                                 ...
                 'filePathOutput','scaling_factor','G_par','par','flatness');
        end
        
% WRITE the clean lfp signal to file
        m = memmapfile(filePathOutput,'Format',{'int16',save_range{2},'data'},'Writable',true);
        m.Data.data(HP, save_range{1}(1,1):save_range{1}(1,2)) = int16(data');
        clear('m');
% WRITE the emg signal to file
        m = memmapfile(filePathEMG,'Format',{'int16',[1 save_range{2}(2)],'data'},'Writable',true);
        tmp_EMG_au = scaling_factor*EMG_au;
        m.Data.data(save_range{1}(1,1):save_range{1}(1,2)) = int16(tmp_EMG_au');
        clear('m');
% $$$     end
end
