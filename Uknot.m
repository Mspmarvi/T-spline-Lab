classdef Uknot<handle
    
    properties
        BackKnot
        NextKnot
        Number
        Value
    end
    
    methods
        function insertAfter(obj,beforeNode)
            obj.BackKnot=beforeNode;
            beforeNode.NextKnot=obj;
        end
        function insertBefore(obj,afterNode)
            obj.NextKnot=afterNode;
            afterNode.BackKnot=obj;
        end
    end
    
end