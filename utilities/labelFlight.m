function Trial = labelFlight(Trial,varargin)
[Stc,overwrite,mode] = DefaultArgs(varargin,{[],false,'get'});

if isempty(Stc),
    Stc = Trial.stc.copy;
end

xyz = Trial.load('xyz');

switch mode
  case 'get'
    sempty = isempty(Stc{'f'});
    if sempty||overwrite,
        if ~sempty
            Stc.states(Stc.gsi('f')) = [];
        end

        load(fullfile(Trial.path.MTAPath,'flight_selection_criteria.mat'));
        fly_fet =  xyz(:,1,3)*usp.pram(:,1)' + repmat(usp.pram(:,2)',xyz.size(1),1)-repmat(xyz(:,end,3),1,size(usp.pram,1));
        
% $$$          figure,hold on,
% $$$          plot(prod(fly_fet<0,2)*100)
% $$$          plot(xyz(:,1,3))
% $$$          plot(xyz(:,end,3),'r')
        
        fper = ThreshCross(prod(fly_fet<0,2),.5,.5*xyz.sampleRate);
% $$$          Lines(fper(:),[],'k');
        Stc.addState(Trial.spath,...
                     Trial.filebase,...
                     fper,...
                     xyz.sampleRate,...
                     Trial.sync.copy,...
                     Trial.sync.data(1),...
                     'flight','f');        
    end

  case 'set'
    figH = figure(48850);
    plot(xyz(:,1,3),xyz(:,end,3),'.');
    xlim(Trial.maze.boundaries(3,:));
    ylim(Trial.maze.boundaries(3,:));
    title('Fit for hyperplane perpendicular to selected dimensions');
    xlabel(['Height of ' xyz.model.Markers{1}.name]);
    ylabel(['Height of ' xyz.model.Markers{end}.name]);
    pram = draw_lines(figH,'line_fit');
    usp.pram   = pram;
    usp.fields = {xyz.model.Markers{1}.name,xyz.model.Markers{end}.name};
    save(fullfile(Trial.path.MTAPath,'flight_selection_criteria.mat'),'usp','-v7.3');

end


Stc.save(1);
Trial.stc = Stc;
Trial.save;
