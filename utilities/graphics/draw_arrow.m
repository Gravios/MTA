function handle = draw_arrow(xcoord,ycoord,varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct("shape",     ...
                 "direction", ...
                 "arrowhead") ...
);
[shape,direction,arrowhead] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

switch shape
  case "curve"
    x = 0;
    y = 0;
    radius = sqrt(sum(xcoord.^2))    
    steps = floor( 2 * pi radius)
    grid = 1:steps)
    angle = grid * 2 * pi / steps
    
    X = x + sin(angle) * radius
    Y = y
  otherwise % line
    points = 
    
    handle = plot(xcoord)
end

    


% UPDATE hash property of Data object
Data.update_hash(DataHash(struct('order',order,'numApplications',numApplications)));



% END MAIN -----------------------------------------------------------------------------------------


% AUX METHODS --------------------------------------------------------------------------------------
% END AUX METHODS ----------------------------------------------------------------------------------


