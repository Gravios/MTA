function [model] = load_patch_model(modelName)

switch modelName
  case 'rat'
    % MODEL - rat
    clear('model');
    % RAT head 
    model.head.patch.color = [0.75,0.75,0.75];
    model.head.patch.vert = ...    
        { ... %   X coords
            [[-0.15],                         ... HeadB
             [-0.30,-0.40,-0.26,-0.18],       ... HeadL
             [-0.11,-0.10, 0.00, 0.10, 0.11], ... Nose
             [ 0.18, 0.26, 0.4, 0.3],         ... HeadR
             [0.15]]                          ... HeadB
            .*50                             ... Scale
            , ... Y coords
            [[-0.15],                         ... HeadB
             [ 0.00, 0.40, 0.90, 1.10],       ... HeadL
             [ 1.20, 1.28, 1.33, 1.28, 1.20], ... Nose
             [ 1.10, 0.90, 0.4, 0.0],         ... HeadR
             [-0.15]]                         ... HeadB
            .*50                             ... Scale
        };

    model.head.midline.color = 'r';
    model.head.midline.vert  = [];
    model.head.overlay.color = [1,0,0];
    model.head.overlay.vert  = [];

    % RAT BODY 
    model.body.patch.color = [0.75,0.75,0.75];
    model.body.patch.vert = ...
        { ... % NeckR         ThoraxR         TTail   ThoraxL              NeckL            Scale   xShift
            [ [ 0.1,  0.60] [ 0.8, 0.85, 0.8],[ 0.0],[-0.8, -0.85, -0.8], [-0.60, -0.1]]  * 48   + 0.0,       ...
            [ [ 0.0, -0.50] [-1.0,-2.00,-2.5],[-2.5],[-2.5, -2.00, -1.0], [-0.50,  0.0]]  * 50   + 0.2        ...
        };
    model.body.midline.color = 'r';
    model.body.midline.vert  = [];
    model.body.overlay.color = [1,0,0];
    model.body.overlay.vert  = [];
  otherwise
    error('MTA:utilities:graphics:load_patch_model:UndefinedModel');
end
