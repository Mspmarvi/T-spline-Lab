classdef ControlPoint<handle
    
    properties
        %
        X;Y;Z;W;
        %
        Uvector;  Vvector;
        %
        RightPoint; LeftPoint; UpPoint; DownPoint;
        %
        RightEdge; LeftEdge; UpEdge; DownEdge
        %
        Number
    end
    properties(Dependent=true)
         CPCU                        %  CentralParameterCoordinateU; %CPCU
         CPCV                         %         CentralParameterCoordinateV; %CPCV
        %
        ParameterSpanU;
        ParameterSpanV;
        %homogenous coordinate
        XW
        YW
        ZW
    end
    
    methods
        function obj=ControlPoint( x, y ,z ,w , uvector , vvector)
            obj.X=x; obj.Y=y;  obj.Z=z; obj.W=w;
            obj.Uvector=uvector;
            obj.Vvector=vvector;
            obj.Number=ControlPoint.counting;
        end
        function cpcv=get.CPCV(obj)
            cpcv=obj.Vvector(3);
        end
        function cpcu=get.CPCU(obj)
            cpcu=obj.Uvector(3);
        end
        function ParameterSpanv=get.ParameterSpanV(obj)
            A=zeros(4,1);
            for i=1:4
                A(i,1)=obj.Vvector(i+1).Value-obj.Vvector(i).Value;
            end
            ParameterSpanv=A;
        end
        function ParameterSpanu=get.ParameterSpanU(obj)
            A=zeros(4,1);
            for i=1:4
                A(i,1)=obj.Uvector(i+1).Value-obj.Uvector(i).Value;
            end
            ParameterSpanu=A;
        end
        function xw=get.XW(obj)
            xw=obj.X*obj.W;
        end
        function yw=get.YW(obj)
            yw=obj.Y*obj.W;
        end
        function zw=get.ZW(obj)
            zw=obj.Z*obj.W;
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