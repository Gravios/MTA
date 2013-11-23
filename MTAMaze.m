classdef MTAMaze
%MTAMaze(name) - contains name and physical parameters of maze
%
    
    properties

        %name - string: 3-6 character name of the maze 
        name

        %shape - string: shape of the maze
        shape

        %visible_volume - numericArray<-double: xyz limits for visualization
        visible_volume

        %boundaries - numericArray<-double: xyz limits for calculations
        boundaries
    end
    
    methods
        function Maze = MTAMaze(name)
            load('MTAMazes.mat');
            for i = 1:length(MTAMazes)
                if strcmp(MTAMazes{i}{1},name),
                    Maze.name = MTAMazes{i}{1};
                    Maze.shape = MTAMazes{i}{2};
                    Maze.visible_volume = MTAMazes{i}{3};
                    Maze.boundaries = MTAMazes{i}{4};
                    return
                elseif i==length(MTAMazes)
                    error('Maze not found in MTAMazes')
                end
            end
        end
    end
            
end

