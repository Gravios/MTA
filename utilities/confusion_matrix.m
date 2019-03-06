function [confmat,labels] =  confusion_matrix(stc1,stc2,varargin)
[labels,sampleRate] = DefaultArgs(varargin,{{},120});

% $$$ if ~isempty(varargin)
% $$$     labels = varargin{1};
% $$$ else 
% $$$     labels = {};
% $$$ end

%% Find common states
s1ep = {};
s2ep = {};
for i = 1:numel(stc1.states),
    for j = 1:numel(stc2.states),
        if strcmp(stc1.states{i}.label,stc2.states{j}.label)...
           &&~isempty(stc1.states{i}.data)...
           &&~isempty(stc2.states{i}.data),
            labels = cat(2,labels,stc1.states{i}.label);
            s1ep{numel(s1ep)+1} = stc1{i,5};
            s2ep{numel(s2ep)+1} = stc2{j,5};
            s1ep{numel(s1ep)}.data(:,2) = s1ep{numel(s1ep)}.data(:,2)-1;
            s2ep{numel(s2ep)}.data(:,2) = s2ep{numel(s2ep)}.data(:,2)-1;
            resample(s1ep{numel(s1ep)}.cast('TimeSeries'),sampleRate);
            resample(s2ep{numel(s2ep)}.cast('TimeSeries'),sampleRate);
            s1ep{numel(s1ep)} = s1ep{numel(s1ep)}.data;
            s2ep{numel(s2ep)} = s2ep{numel(s2ep)}.data;
        end
    end
end

%assert(~isempty(s1ind)&~isempty(s2ind),'No matching States')



confmat = confusionmat(sum(cell2mat(s1ep).*repmat(1:numel(labels),[size(s1ep{1},1),1]),2),...
                       sum(cell2mat(s2ep).*repmat(1:numel(labels),[size(s2ep{1},1),1]),2));
