classdef Edge<handle
  
    
    properties
        Direction
        FirstTip
        LastTip
        FirstControlPoint
        LastControlPoint
        DirectionParameterCoordinate
        Number
    end
    
    methods
        function obj=Edge(direction)
            obj.Direction=direction;
            obj.Number=Edge.counting;
        end
    end
    methods(Static=true)
        function number=counting
            persistent i
            if isempty(i)
                i=1;
            else
                i=i+1;
            end
            number=i;
        end
    end
    
    
end