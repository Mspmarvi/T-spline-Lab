classdef T_Surface<handle
    
    properties%(all of them are array of objects)
        controlpoint %vector of ControlPoint objects
        edge %vector of Edge objects
        U   %vector of Uknot objects
        V   %vector of Vknot objects
        C0_points   
        C0_curveU  
        C0_curveV  
    end
    
    methods
        function obj=T_Surface(varargin)%we have two types of input, type 1: no inputs that means you construct a T-surface with no details, 
            %type 2:  matrix by size of (n+1,m+1,4)
            %data is a matrix by size of (n+1,m+1,4); dimension 3 indicates x, y, z, w;
            %which represents a grid of NURBS controlpoints
            if nargin==1
                data=varargin{1};
                %% controlpoint construction
                nplus=size(data,1);
                n=nplus-1;
                mplus=size(data,2);
                m=mplus-1;
                %clamped boundry condition to construct the u&vVector
                degree=3; %in T-Spline,degree is always 3
                U=zeros((n+degree+1)+1,1);
                V=zeros((m+degree+1)+1,1);
                %%
                %clamped condition for knot vectors
%                 du=1/(n+degree-2*degree+1);
%                 dv=1/(m+degree-2*degree+1);
%                 for i=degree+1:n+1
%                     U(i+1)=U(i)+du;
%                 end
%                 for i=degree+1:m+1
%                     V(i+1)=V(i)+dv;
%                 end
%                 U(n+2:end,1)=1;
%                 V(m+2:end,1)=1;
                %
% uniform condition for knot vectors 
%%
                du=1/(n+degree+1);
                dv=1/(m+degree+1);

                for i=1:n+degree+1
                    U(i+1)=U(i)+du;
                end
                U(end) = 1;
                for i=1:m+degree+1
                    V(i+1)=V(i)+dv;
                end
                V(end) = 1;
                  %%
                  
                  
                  
                obj.U=Uknot.empty(n+degree+1+1,0);
                obj.V=Vknot.empty(m+degree+1+1,0);
                
                %allocating the U & V
                for i=1:n+degree+1+1
                    obj.U(i).Value=U(i);
                    obj.U(i).Number=i;
                    if i>1
                        obj.U(i).BackKnot=obj.U(i-1);
                        obj.U(i-1).NextKnot=obj.U(i);
                    end
                    
                end
                
                for i=1:m+degree+1+1
                    obj.V(i).Value=V(i);
                    obj.V(i).Number=i;
                    if i>1
                        obj.V(i).BackKnot=obj.V(i-1);
                        obj.V(i-1).NextKnot=obj.V(i);
                    end
                    
                end
                %%
                %put controlpoint objects in cell to create default mesh a
                %NURBS surface
                cp=cell(n+1,m+1);
                c=0;
                obj.controlpoint=ControlPoint.empty((n+1)*(m+1),0);
                for i=1:n+1
                    for j=1:m+1
                        c=c+1;
                        cp{i,j}=ControlPoint(data(i,j,1),data(i,j,2),data(i,j,3),data(i,j,4),obj.U(i:i+degree+1),obj.V(j:j+degree+1));
                        obj.controlpoint(c)=cp{i,j};
                    end
                end
                
                %neighberhood(put neighbering controlpoints to controlpoint object)
                for i=1:n+1
                    for j=1:m+1
                        if i>1
                            cp{i,j}.UpPoint=cp{i-1,j};
                        end
                        if i<n+1
                            cp{i,j}.DownPoint=cp{i+1,j};
                        end
                        if j>1
                            cp{i,j}.LeftPoint=cp{i,j-1};
                        end
                        if j<m+1
                            cp{i,j}.RightPoint=cp{i,j+1};
                        end
                    end
                end
                
                %% edge construction in cells
                edH=cell(n+1,m+2);
                edV=cell(n+2,m+1);
                c=0;
                %Horizontal
                obj.edge=Edge.empty((n+1)*(m+2)+(n+2)*(m+1),0);
                for i=1:n+1
                    for j=1:m+2
                        edH{i,j}=Edge('horizontal');
                        c=c+1;
                        obj.edge(c)=edH{i,j};
                        if j>1
                            edH{i,j}.FirstControlPoint=cp{i,j-1};
                        end
                        if j<m+2
                            edH{i,j}.LastControlPoint=cp{i,j};
                        end
                        edH{i,j}.DirectionParameterCoordinate=obj.U(i+2);
                        edH{i,j}.FirstTip=obj.V(j+1);
                        edH{i,j}.LastTip=obj.V(j+2);
                    end
                end
                %Vertical
                for i=1:n+2
                    for j=1:m+1
                        c=c+1;
                        edV{i,j}=Edge('vertical');
                        obj.edge(c)=edV{i,j};
                        if i>1
                            edV{i,j}.FirstControlPoint=cp{i-1,j};
                        end
                        if i<n+2
                            edV{i,j}.LastControlPoint=cp{i,j};
                        end
                        edV{i,j}.DirectionParameterCoordinate=obj.V(j+2);
                        edV{i,j}.FirstTip=obj.U(i+1);
                        edV{i,j}.LastTip=obj.U(i+2);
                    end
                end
                
                
                %neighberhood(put neighbering controlpoints to controlpoint object)
                for i=1:n+1
                    for j=1:m+1
                        cp{i,j}.RightEdge=edH{i,j+1};
                        cp{i,j}.LeftEdge=edH{i,j};
                        cp{i,j}.UpEdge=edV{i,j};
                        cp{i,j}.DownEdge=edV{i+1,j};
                    end
                end
            else
                %%this is used when you don't have inputs for surface
                obj.controlpoint=[];
                obj.U=[];
                obj.V=[];
                obj.edge=[];
            end
        end
        %%
        function plotmesh(obj) %this fucntion is to plot control grid in cartesian. the input is just T-spline surface 
            figure(1)
            hold on
            %%
            try
                %this part is to write a the number of edges in plot
                for i=1:numel(obj.edge)
                    if ~isempty(obj.edge(i).FirstControlPoint)  &&  ~isempty(obj.edge(i).LastControlPoint)
%                         text((obj.edge(i).FirstControlPoint.X+obj.edge(i).LastControlPoint.X)/2,....
%                             (obj.edge(i).FirstControlPoint.Y+obj.edge(i).LastControlPoint.Y)/2,...
%                             (obj.edge(i).FirstControlPoint.Z+obj.edge(i).LastControlPoint.Z)/2,....
%                             num2str(obj.edge(i).Number),'Color',[0 0 1],'Fontsize',11);
                    end
                end
            catch
                error(['i= ',num2str(i)]);
            end
            
            try
                %this part is to draw edges in plot
                for i=1:numel(obj.edge)
                    if ~isempty(obj.edge(i).FirstControlPoint)  &&  ~isempty(obj.edge(i).LastControlPoint)
                        line([obj.edge(i).FirstControlPoint.X obj.edge(i).LastControlPoint.X],....
                            [obj.edge(i).FirstControlPoint.Y obj.edge(i).LastControlPoint.Y],....
                            [obj.edge(i).FirstControlPoint.Z obj.edge(i).LastControlPoint.Z],'linewidth',2,'color',[0.5 0 1]);
                    end
                end
            catch
                error(['i= ',num2str(i)]);
            end
            %%
            try
                %this part is draw control points in plot
                for i=1:numel(obj.controlpoint)
                    plot3(obj.controlpoint(i).X,obj.controlpoint(i).Y,obj.controlpoint(i).Z,'o','markersize',3,'markerfacecolor','k','markeredgecolor','g');
                end
            catch
                error(['i= ',num2str(i)]);
            end
            %%
            try
                %this part is to draw control point in plot
                for i=1:numel(obj.controlpoint)
%                     text(obj.controlpoint(i).X,obj.controlpoint(i).Y,obj.controlpoint(i).Z,num2str(obj.controlpoint(i).Number),'Fontsize',12,'color','k')
                end
            catch
                error(['i= ',num2str(i)]);
            end
            axis equal
            view(20,20);
        end
        %%
        function plot(obj) % this part is to plot surface point by point.the input is T-spline surface
            figure(1)
            hold on
            p=3;% degree
            s=0 ; %save the coordinate in to matrix
            for u=obj.U(p+1).Value:0.01:obj.U(end-p).Value %it must be 0.005
                for v=obj.V(p+1).Value:0.02:obj.V(end-p).Value %it must be 0.005
                    
                    for i=1:numel(obj.controlpoint)
                        if (u>=obj.controlpoint(i).Uvector(1).Value && u<obj.controlpoint(i).Uvector(end).Value)
                            if v>=obj.controlpoint(i).Vvector(1).Value && v<obj.controlpoint(i).Vvector(end).Value
                                
                                [N,Q]=T_Surface.basisfunction(obj.controlpoint(i),u,v,p);
                                X=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).X*obj.controlpoint(i).W;
                                Y=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).Y*obj.controlpoint(i).W;
                                Z=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).Z*obj.controlpoint(i).W;
                                Weight=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).W;
                                if ~exist('x','var')
                                    x=X;
                                    y=Y;
                                    z=Z;
                                    w=Weight;
                                else
                                    x=x+X;
                                    y=y+Y;
                                    z=z+Z;
                                    w=w+Weight;
                                end
                            end
                        end
                    end
                    hold on
                    if exist('x','var')
                        %%
                        s=s+1;
                        coordinateX(s)=x/w;
                        coordinateY(s)=y/w;
                        coordinateZ(s)=z/w;
                        %%
                        clear x y z
                    end
                end
            end
            plot3(coordinateX(:),coordinateY(:),coordinateZ(:),'mo','markerfacecolor','g','markersize',2,'markeredgecolor','g')
            axis equal
            view(20,20);
        end
        %%
        function plotSurf(obj) %this part is to plot surface. the input is T-spline surface
            figure(1)
            hold on
            data = extract_Data(obj);
            surf(data(:,:,1),data(:,:,2),data(:,:,3));
        end
        %%
        function obj=InsertPoint(obj,edgenumber,varargin) % this function is used to add control points by the algorithm which is introduced in 2003 by sederberg
            %%
            %there are three types of using this function
                                %% FIRST %%
            % obj = InsertPoint(obj,edgenumber)  
            % to simplified the using of inserting control point user can
            % add T-spline object and the number of edge that user wants to
            % split then add a control point. parameter of new control point
            % is in the middle of former(selected) edge
                                %% SECOND %%
             % obj = InsertPoint(obj,edgenumber,knot)
             % the inputs are T-spline object, the number of edge that user
             % wants to split and a existed knot which it's parameter is
             % locatad on edge.
             % this type of inputs is mostly used by program and functions in some
             % situations
                                %% THIRD %%
             % obj = InsertPoint(obj,edgenumber,value)
             % the inputs are like previous type except value which is the
             % parameter value
             %this type mostly used by other programs and functions
             %%
             
             
            edgepopulation=numel(obj.edge);
            %to find the edge(the variable 'i' is important that indicate
            %the index of edge number(as amature programmer it is programmed : )  )
            for i=1:edgepopulation
                if obj.edge(i).Number==edgenumber
                    SelectedEdge=obj.edge(i); %SelectedEdge is edge you choose to split
                    break
                end
            end
            %
            edgeOrientparameter=SelectedEdge.DirectionParameterCoordinate;    %edgeOrientparameter is the parameter of along woth of SelectedEdge  
            direct=SelectedEdge.Direction; %direct is direction of SelectedEdge ('horizontal' or 'vertical')
            
                                %% newknot %%
            %this is to make new knot with parameter we add as input and add it to U or V vector of T-spline object 
            if nargin==2
                switch direct
                    case 'horizontal'
                        if SelectedEdge.LastTip.Number-SelectedEdge.FirstTip.Number==2
                            number=SelectedEdge.FirstTip.Number+1;
                            newknot=findobj(obj.V,'Number',number);
                        else
                            newknot=Vknot;
                            insertAfter(newknot,SelectedEdge.FirstTip);
                            insertBefore(newknot,SelectedEdge.LastTip);
                            %resorting in V
                            c=1;
                            knot=obj.V(1);
                            vv=Vknot.empty(numel(obj.V)+1,0);
                            vv(1)=knot;
                            while(true)
                                knot=knot.NextKnot;
                                if c<=numel(obj.V)
                                    vv(c+1)=knot;
                                    knot.Number=c+1;
                                    c=c+1;
                                else
                                    break
                                end
                            end
                            obj.V=vv;
                        end
                    case 'vertical'
                        if SelectedEdge.LastTip.Number-SelectedEdge.FirstTip.Number==2
                            number=SelectedEdge.FirstTip.Number+1;
                            newknot=findobj(obj.U,'Number',number);
                        else
                            newknot=Uknot;
                            insertAfter(newknot,SelectedEdge.FirstTip);
                            insertBefore(newknot,SelectedEdge.LastTip);
                            %resorting in U
                            c=1;
                            knot=obj.U(1);
                            uu=Uknot.empty(numel(obj.U)+1,0);
                            uu(1)=knot;
                            while(true)
                                knot=knot.NextKnot;
                                if c<=numel(obj.U)
                                    uu(c+1)=knot;
                                    knot.Number=c+1;
                                    c=c+1;
                                else
                                    break
                                end
                            end
                            obj.U=uu;
                        end
                end
                newknot.Value=(SelectedEdge.FirstTip.Value +SelectedEdge.LastTip.Value )/2;
            elseif nargin==3
                newknot=varargin{1};%%IT MIGHT BE WRONG
            elseif nargin==4
                switch direct
                    case 'horizontal'
                        newknot=Vknot;
                        k=SelectedEdge.FirstTip.Number;
                        j=SelectedEdge.LastTip.Number;
                        value=varargin{2};
                        while(true)
                            if value<obj.V(k+1).Value && value>=obj.V(k).Value
                                newknot.Number=k+1;
                                newknot.Value=value;
                                v=Vknot.empty(numel(obj.V)+1,0);
                                v(1:k)=obj.V(1:k);
                                v(k+1)=newknot;
                                v(k+2:numel(obj.V)+1)=obj.V(k+1:end);
                                obj.V=v;
                                break
                            end
                            k=k+1;
                        end
                        insertAfter(newknot,obj.V(k));
                        insertBefore(newknot,obj.V(k+1));
                        %resorting in V
                        for c=1:numel(obj.V)
                            obj.V(c).Number=c;
                        end
                        
                    case 'vertical'
                        newknot=Uknot;
                        k=SelectedEdge.FirstTip.Number;
                        j=SelectedEdge.LastTip.Number;
                        value=varargin{2};
                        while(true)
                            if value<obj.U(k+1).Value && value>=obj.U(k).Value
                                newknot.Number=k+1;
                                newknot.Value=value;
                                u=Uknot.empty(numel(obj.U)+1,0);
                                u(1:k)=obj.U(1:k);
                                u(k+1)=newknot;
                                u(k+2:numel(obj.U)+1)=obj.U(k+1:end);
                                obj.U=u;
                                break
                            end
                            k=k+1;
                        end
                        insertAfter(newknot,obj.U(k));
                        insertBefore(newknot,obj.U(k+1));
                        %resorting in U
                        for c=1:numel(obj.U)
                            obj.U(c).Number=c;
                        end
                end
            end
            %%
                                                    %%  neighboring edges of Selected edge %%
            % this part is to check neighboring edges of Selected edge to
            % find out that if they exist or not
            % if they don't exist we most build a edge on every side of
            % Selected edge in same orientation of Selected edge
            switch direct
                
                case  'horizontal'
                    if isempty(SelectedEdge.FirstControlPoint.LeftEdge) % to find if there is a edge on left of the Selected edge. if there is no left edge we most build it
                        situation=false;
                        
                        for c=SelectedEdge.FirstTip.Number-1:-1:1
                            group=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(c));
                            for j=1:numel(group)
                                if group(j).FirstTip.Number<SelectedEdge.DirectionParameterCoordinate.Number && .....
                                        group(j).LastTip.Number>SelectedEdge.DirectionParameterCoordinate.Number
                                    situation=true;
                                    obj.InsertPoint(group(j).Number,SelectedEdge.DirectionParameterCoordinate);
                                end
                            end
                            if situation
                                break
                            end
                        end
                    end
                    if isempty(SelectedEdge.LastControlPoint.RightEdge) % to find if there is a edge on right of the Selected edge. if there is no left edge we most to build it
                        situation=false;
                        for c=SelectedEdge.LastTip.Number+1:numel(obj.V)
                            group=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(c));
                            for j=1:numel(group)
                                if group(j).FirstTip.Number<SelectedEdge.DirectionParameterCoordinate.Number && .....
                                        group(j).LastTip.Number>SelectedEdge.DirectionParameterCoordinate.Number
                                    situation=true;
                                    obj.InsertPoint(group(j).Number,SelectedEdge.DirectionParameterCoordinate);
                                end
                            end
                            if situation
                                break
                            end
                        end
                    end
                case  'vertical' % it's the same approach in the previous direction
                    if isempty(SelectedEdge.FirstControlPoint.UpEdge)
                        situation=false;
                        for c=SelectedEdge.FirstTip.Number-1:-1:1
                            group=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(c));
                            for j=1:numel(group)
                                if group(j).FirstTip.Number<SelectedEdge.DirectionParameterCoordinate.Number && .....
                                        group(j).LastTip.Number>SelectedEdge.DirectionParameterCoordinate.Number
                                    situation=true;
                                    obj.InsertPoint(group(j).Number,SelectedEdge.DirectionParameterCoordinate);
                                end
                            end
                            if situation
                                break
                            end
                        end
                    end
                    if isempty(SelectedEdge.LastControlPoint.DownEdge)
                        situation=false;
                        for c=SelectedEdge.LastTip.Number+1:numel(obj.U)
                            group=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(c));
                            for j=1:numel(group)
                                if group(j).FirstTip.Number<SelectedEdge.DirectionParameterCoordinate.Number && .....
                                        group(j).LastTip.Number>SelectedEdge.DirectionParameterCoordinate.Number
                                    situation=true;
                                    obj.InsertPoint(group(j).Number,SelectedEdge.DirectionParameterCoordinate);
                                end
                            end
                            if situation
                                break
                            end
                        end
                    end
            end
            
            
                                    %% domain checking %%
            % this part is to check domian of neighboring control points 
            %if there are some inequality of domains we most add control
            %points to satisfy the pre-condition of adding control point on Selected edge 
            switch direct
                case 'vertical'
                    %one step down domain checking
                    for localknotnumber=2:-1:1
                        
                        for iteration=1:3
                            %%
                            if ~isempty(SelectedEdge.FirstControlPoint.UpPoint)
                                if SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber).Number<SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.LeftEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.LeftPoint.LeftEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Vvector(localknotnumber))
                                elseif  SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber).Number>obj.edge(i).FirstControlPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.LeftEdge.Number;
                                    end
                                    obj.InsertPoint( number_of_edge,SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber))
                                end
                            end
                            %%
                            if SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number<SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number
                                if localknotnumber==2
                                    number_of_edge=SelectedEdge.FirstControlPoint.LeftEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.LeftEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Vvector(localknotnumber))
                            elseif SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number>SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number
                                if localknotnumber==2
                                    number_of_edge=SelectedEdge.LastControlPoint.LeftEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.LastControlPoint.LeftPoint.LeftEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Vvector(localknotnumber))
                            end
                            %%
                            if ~isempty(SelectedEdge.LastControlPoint.DownPoint)
                                if SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number<SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.LastControlPoint.LeftEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.LeftPoint.LeftEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber))
                                elseif SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number>SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.LastControlPoint.DownPoint.LeftEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.DownPoint.LeftPoint.LeftEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Vvector(localknotnumber))
                                end
                            end
                        end
                        
                    end
                    for localknotnumber=4:1:5
                        
                        for iteration=1:3
                            %%
                            if ~isempty(SelectedEdge.FirstControlPoint.UpPoint)
                                if SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber).Number>SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.RightEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.RightPoint.RightEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Vvector(localknotnumber))
                                elseif  SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber).Number<SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.FirstControlPoint.RightEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.RightPoint.RightEdge.Number;
                                    end
                                    obj.InsertPoint( number_of_edge,SelectedEdge.FirstControlPoint.UpPoint.Vvector(localknotnumber))
                                end
                            end
                            %%
                            if SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number>SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number
                                if localknotnumber==4
                                    number_of_edge=SelectedEdge.FirstControlPoint.RightEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.FirstControlPoint.RightPoint.RightEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Vvector(localknotnumber))
                            elseif SelectedEdge.FirstControlPoint.Vvector(localknotnumber).Number<SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number
                                if localknotnumber==4
                                    number_of_edge=SelectedEdge.LastControlPoint.RightEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.LastControlPoint.RightPoint.RightEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Vvector(localknotnumber))
                            end
                            %%
                            if ~isempty(SelectedEdge.LastControlPoint.DownPoint)
                                if SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number>SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.LastControlPoint.RightEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.RightPoint.RightEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber))
                                elseif SelectedEdge.LastControlPoint.Vvector(localknotnumber).Number<SelectedEdge.LastControlPoint.DownPoint.Vvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.LastControlPoint.DownPoint.RightEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.DownPoint.RightPoint.RightEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Vvector(localknotnumber))
                                end
                            end
                        end
                        
                    end
                case 'horizontal'
                    for localknotnumber=2:-1:1
                        
                        for iteration=1:3
                            %%
                            if ~isempty(SelectedEdge.FirstControlPoint.LeftPoint)
                                if SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber).Number < SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.UpEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.UpPoint.UpEdge.Number;
                                    end
                                    
                                    obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Uvector(localknotnumber))
                                elseif  SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber).Number>SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.UpEdge.Number;
                                    end
                                    obj.InsertPoint( number_of_edge,SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber))
                                end
                            end
                            %%
                            if SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number<SelectedEdge.LastControlPoint.Uvector(localknotnumber).Number
                                if localknotnumber==2
                                    number_of_edge=SelectedEdge.FirstControlPoint.UpEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.FirstControlPoint.UpPoint.UpEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Uvector(localknotnumber))
                            elseif SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number>SelectedEdge.LastControlPoint.Uvector(localknotnumber).Number
                                if localknotnumber==2
                                    number_of_edge=SelectedEdge.LastControlPoint.UpEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.LastControlPoint.UpPoint.UpEdge.Number;
                                end
                                
                                obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Uvector(localknotnumber))
                            end
                            %%
                            if ~isempty(SelectedEdge.LastControlPoint.RightPoint)
                                if SelectedEdge.LastControlPoint.Uvector(localknotnumber).Number<SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.LastControlPoint.UpEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.UpPoint.UpEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber))
                                elseif SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number>SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==2
                                        number_of_edge=SelectedEdge.LastControlPoint.RightPoint.UpEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.RightPoint.UpPoint.UpEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Uvector(localknotnumber))
                                end
                            end
                        end
                        
                    end
                    for localknotnumber=4:1:5
                        
                        for iteration=1:3
                            %%
                            if ~isempty(SelectedEdge.FirstControlPoint.LeftPoint)
                                if SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber).Number>SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.DownEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.LeftPoint.DownPoint.DownEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,obj.edge(i).FirstControlPoint.Uvector(localknotnumber))
                                elseif  SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber).Number<SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.FirstControlPoint.DownEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.FirstControlPoint.DownPoint.DownEdge.Number;
                                    end
                                    obj.InsertPoint( number_of_edge,SelectedEdge.FirstControlPoint.LeftPoint.Uvector(localknotnumber))
                                end
                            end
                            %%
                            if SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number>obj.edge(i).LastControlPoint.Uvector(localknotnumber).Number
                                if localknotnumber==4
                                    number_of_edge=SelectedEdge.FirstControlPoint.DownEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.FirstControlPoint.DownPoint.DownEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Uvector(localknotnumber))
                            elseif SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number<SelectedEdge.LastControlPoint.Uvector(localknotnumber).Number
                                if localknotnumber==4
                                    number_of_edge=SelectedEdge.LastControlPoint.DownEdge.Number;
                                else
                                    number_of_edge=SelectedEdge.LastControlPoint.DownPoint.DownEdge.Number;
                                end
                                obj.InsertPoint(number_of_edge,SelectedEdge.FirstControlPoint.Uvector(localknotnumber))
                            end
                            %%
                            if ~isempty(SelectedEdge.LastControlPoint.RightPoint)
                                if SelectedEdge.LastControlPoint.Uvector(localknotnumber).Number>SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.LastControlPoint.DownEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.DownPoint.DownEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber))
                                elseif SelectedEdge.FirstControlPoint.Uvector(localknotnumber).Number<SelectedEdge.LastControlPoint.RightPoint.Uvector(localknotnumber).Number
                                    if localknotnumber==4
                                        number_of_edge=SelectedEdge.LastControlPoint.RightPoint.DownEdge.Number;
                                    else
                                        number_of_edge=SelectedEdge.LastControlPoint.RightPoint.DownPoint.DownEdge.Number;
                                    end
                                    obj.InsertPoint(number_of_edge,SelectedEdge.LastControlPoint.Uvector(localknotnumber))
                                end
                            end
                        end
                        
                    end
            end
            
            %create new controlpoint (according to sederberg's paper)
            switch direct
                case 'horizontal'
                    P1= SelectedEdge.FirstControlPoint.LeftPoint;
                    P2=SelectedEdge.FirstControlPoint;
                    P4=SelectedEdge.LastControlPoint;
                    P5=SelectedEdge.LastControlPoint.RightPoint;
                    
                    %homogenous coordinates
                    P2XW=(P1.XW*(P4.CPCV.Value-newknot.Value)+(newknot.Value-P1.Vvector(2).Value)*P2.XW)/(P4.CPCV.Value-P1.Vvector(2).Value);
                    P2YW=(P1.YW*(P4.CPCV.Value-newknot.Value)+(newknot.Value-P1.Vvector(2).Value)*P2.YW)/(P4.CPCV.Value-P1.Vvector(2).Value);
                    P2ZW=(P1.ZW*(P4.CPCV.Value-newknot.Value)+(newknot.Value-P1.Vvector(2).Value)*P2.ZW)/(P4.CPCV.Value-P1.Vvector(2).Value);
                    P2W=(P1.W*(P4.CPCV.Value-newknot.Value)+(newknot.Value-P1.Vvector(2).Value)*P2.W)/(P4.CPCV.Value-P1.Vvector(2).Value);
                    
                    P4XW=(P4.XW*(P5.Vvector(4).Value-newknot.Value)+P5.XW*(newknot.Value-P2.CPCV.Value))/(P5.Vvector(4).Value-P2.CPCV.Value);
                    P4YW=(P4.YW*(P5.Vvector(4).Value-newknot.Value)+P5.YW*(newknot.Value-P2.CPCV.Value))/(P5.Vvector(4).Value-P2.CPCV.Value);
                    P4ZW=(P4.ZW*(P5.Vvector(4).Value-newknot.Value)+P5.ZW*(newknot.Value-P2.CPCV.Value))/(P5.Vvector(4).Value-P2.CPCV.Value);
                    P4W=(P4.W*(P5.Vvector(4).Value-newknot.Value)+P5.W*(newknot.Value-P2.CPCV.Value))/(P5.Vvector(4).Value-P2.CPCV.Value);
                    
                    P3XW=(P2.XW*(P5.CPCV.Value-newknot.Value)+P4.XW*(newknot.Value-P1.CPCV.Value))/(P5.CPCV.Value-P1.CPCV.Value);
                    P3YW=(P2.YW*(P5.CPCV.Value-newknot.Value)+P4.YW*(newknot.Value-P1.CPCV.Value))/(P5.CPCV.Value-P1.CPCV.Value);
                    P3ZW=(P2.ZW*(P5.CPCV.Value-newknot.Value)+P4.ZW*(newknot.Value-P1.CPCV.Value))/(P5.CPCV.Value-P1.CPCV.Value);
                    P3W=(P2.W*(P5.CPCV.Value-newknot.Value)+P4.W*(newknot.Value-P1.CPCV.Value))/(P5.CPCV.Value-P1.CPCV.Value);
                    
                    %Non-homogeneous coordinate
                    P2.X=P2XW/P2W;
                    P2.Y=P2YW/P2W;
                    P2.Z=P2ZW/P2W;
                    P2.W=1;
                    
                    P3X=P3XW/P3W;
                    P3Y=P3YW/P3W;
                    P3Z=P3ZW/P3W;
                    P3W=1;
                    
                    P4.X=P4XW/P4W;
                    P4.Y=P4YW/P4W;
                    P4.Z=P4ZW/P4W;
                    P4.W=1;
                case 'vertical'
                    P1= SelectedEdge.FirstControlPoint.UpPoint;
                    P2=SelectedEdge.FirstControlPoint;
                    P4=SelectedEdge.LastControlPoint;
                    P5=SelectedEdge.LastControlPoint.DownPoint;
                    
                    %homogenous coordinates
                    P2XW=(P1.XW*(P4.CPCU.Value-newknot.Value)+(newknot.Value-P1.Uvector(2).Value)*P2.XW)/(P4.CPCU.Value-P1.Uvector(2).Value);
                    P2YW=(P1.YW*(P4.CPCU.Value-newknot.Value)+(newknot.Value-P1.Uvector(2).Value)*P2.YW)/(P4.CPCU.Value-P1.Uvector(2).Value);
                    P2ZW=(P1.ZW*(P4.CPCU.Value-newknot.Value)+(newknot.Value-P1.Uvector(2).Value)*P2.ZW)/(P4.CPCU.Value-P1.Uvector(2).Value);
                    P2W=(P1.W*(P4.CPCU.Value-newknot.Value)+(newknot.Value-P1.Uvector(2).Value)*P2.W)/(P4.CPCU.Value-P1.Uvector(2).Value);
                    
                    P4XW=(P4.XW*(P5.Uvector(4).Value-newknot.Value)+P5.XW*(newknot.Value-P2.CPCU.Value))/(P5.Uvector(4).Value-P2.CPCU.Value);
                    P4YW=(P4.YW*(P5.Uvector(4).Value-newknot.Value)+P5.YW*(newknot.Value-P2.CPCU.Value))/(P5.Uvector(4).Value-P2.CPCU.Value);
                    P4ZW=(P4.ZW*(P5.Uvector(4).Value-newknot.Value)+P5.ZW*(newknot.Value-P2.CPCU.Value))/(P5.Uvector(4).Value-P2.CPCU.Value);
                    P4W=(P4.W*(P5.Uvector(4).Value-newknot.Value)+P5.W*(newknot.Value-P2.CPCU.Value))/(P5.Uvector(4).Value-P2.CPCU.Value);
                    
                    P3XW=(P2.XW*(P5.CPCU.Value-newknot.Value)+P4.XW*(newknot.Value-P1.CPCU.Value))/(P5.CPCU.Value-P1.CPCU.Value);
                    P3YW=(P2.YW*(P5.CPCU.Value-newknot.Value)+P4.YW*(newknot.Value-P1.CPCU.Value))/(P5.CPCU.Value-P1.CPCU.Value);
                    P3ZW=(P2.ZW*(P5.CPCU.Value-newknot.Value)+P4.ZW*(newknot.Value-P1.CPCU.Value))/(P5.CPCU.Value-P1.CPCU.Value);
                    P3W=(P2.W*(P5.CPCU.Value-newknot.Value)+P4.W*(newknot.Value-P1.CPCU.Value))/(P5.CPCU.Value-P1.CPCU.Value);
                    
                    %Non-homogeneous coordinate
                    P2.X=P2XW/P2W;
                    P2.Y=P2YW/P2W;
                    P2.Z=P2ZW/P2W;
                    P2.W=1;
                    
                    P3X=P3XW/P3W;
                    P3Y=P3YW/P3W;
                    P3Z=P3ZW/P3W;
                    P3W=1;
                    
                    P4.X=P4XW/P4W;
                    P4.Y=P4YW/P4W;
                    P4.Z=P4ZW/P4W;
                    P4.W=1;
            end
            
            %%
            % create new edge object
            %we have 2 new edge after splitting. the first piece is created
            %by modifying the older edge and the second one created by creating new edge object
            %newedge is second edge
            %in this section also we will create a new controlpoint that we
            %intented to add at the first place
            
            newEdge=Edge(direct);
            newEdge.FirstTip=newknot;
            newEdge.LastTip=SelectedEdge.LastTip;
            
            switch direct
                case 'horizontal'
                    newControlPoint=ControlPoint(P3X,P3Y,P3Z,P3W,SelectedEdge.FirstControlPoint.Uvector,[SelectedEdge.FirstControlPoint.LeftPoint.CPCV;....
                        SelectedEdge.FirstControlPoint.CPCV;newknot;SelectedEdge.LastControlPoint.CPCV;SelectedEdge.LastControlPoint.RightPoint.CPCV]);
    
                    %other ControlPoints
                    %it start first point from the left
                    SelectedEdge.FirstControlPoint.LeftPoint.Vvector(5)=newknot;
                    
                    SelectedEdge.FirstControlPoint.RightPoint=newControlPoint;
                    SelectedEdge.FirstControlPoint.Vvector(4)=newknot;
                    SelectedEdge.FirstControlPoint.Vvector(5)=SelectedEdge.LastControlPoint.CPCV;
                    
                    SelectedEdge.LastControlPoint.LeftPoint=newControlPoint;
                    SelectedEdge.LastControlPoint.Vvector(1)=SelectedEdge.FirstControlPoint.CPCV;
                    SelectedEdge.LastControlPoint.Vvector(2)=newknot;
                    
                    SelectedEdge.LastControlPoint.RightPoint.Vvector(1)=newknot;
                    
                    %modifying new edge
                    
                    newEdge.FirstControlPoint=newControlPoint;
                    newEdge.LastControlPoint=SelectedEdge.LastControlPoint;
                    newEdge.DirectionParameterCoordinate=SelectedEdge.DirectionParameterCoordinate;
                    %modifying older edge
                    SelectedEdge.LastTip=newknot;
                    SelectedEdge.LastControlPoint=newControlPoint;
                    
                    %resorting the edge
                    n=numel(obj.edge);
                    edge_array=Edge.empty(n+1,0);
                    edge_array=obj.edge(:);
                    
                    obj.edge=edge_array;
                    obj.edge(n+1)=newEdge;
                    %resorting the controlpoint
                    
                    n=max(size(obj.controlpoint));
                    controlpoint_array=ControlPoint.empty(n+1,0);
                    controlpoint_array=obj.controlpoint(:);
                    controlpoint_array(n+1)=newControlPoint;
                    
                    obj.controlpoint=controlpoint_array;
                    
                    newControlPoint.LeftPoint=SelectedEdge.FirstControlPoint;
                    newControlPoint.RightPoint=newEdge.LastControlPoint;
                    newControlPoint.LeftEdge=SelectedEdge;
                    newControlPoint.RightEdge=newEdge;
                case 'vertical'
                    newControlPoint=ControlPoint(P3X,P3Y,P3Z,P3W,[SelectedEdge.FirstControlPoint.UpPoint.CPCU;....
                        SelectedEdge.FirstControlPoint.CPCU;newknot;SelectedEdge.LastControlPoint.CPCU;SelectedEdge.LastControlPoint.DownPoint.CPCU]....
                        ,SelectedEdge.FirstControlPoint.Vvector);
                    %other ControlPoints
                    %it start first point from the left
                    SelectedEdge.FirstControlPoint.UpPoint.Uvector(5)=newknot;
                    
                    SelectedEdge.FirstControlPoint.DownPoint=newControlPoint;
                    SelectedEdge.FirstControlPoint.Uvector(4)=newknot;
                    SelectedEdge.FirstControlPoint.Uvector(5)=SelectedEdge.LastControlPoint.CPCU;
                    
                    SelectedEdge.LastControlPoint.UpPoint=newControlPoint;
                    SelectedEdge.LastControlPoint.Uvector(1)=SelectedEdge.FirstControlPoint.CPCU;
                    SelectedEdge.LastControlPoint.Uvector(2)=newknot;
                    
                    SelectedEdge.LastControlPoint.DownPoint.Uvector(1)=newknot;
                    
                    %modifying new edge
                    
                    newEdge.FirstControlPoint=newControlPoint;
                    newEdge.LastControlPoint=SelectedEdge.LastControlPoint;
                    newEdge.DirectionParameterCoordinate=SelectedEdge.DirectionParameterCoordinate;
                    %modifying older edge
                    SelectedEdge.LastTip=newknot;
                    SelectedEdge.LastControlPoint=newControlPoint;
                    
                    %resorting the edge
                    n=numel(obj.edge);
                    edge_array=Edge.empty(n+1,0);
                    edge_array=obj.edge(:);
                    
                    obj.edge=edge_array;
                    obj.edge(n+1)=newEdge;
                    
                    %resorting the controlpoint
                    
                    n=numel(obj.controlpoint);
                    controlpoint_array=ControlPoint.empty(n+1,0);
                    controlpoint_array(1:n)=obj.controlpoint(1:n);
                    controlpoint_array(n+1)=newControlPoint;
                    obj.controlpoint=controlpoint_array;
                    
                    newControlPoint.UpPoint=SelectedEdge.FirstControlPoint;
                    newControlPoint.DownPoint=newEdge.LastControlPoint;
                    newControlPoint.UpEdge=SelectedEdge;
                    newControlPoint.DownEdge=newEdge;
                    
            end
            
                                %% T-spline rules %%
            %this section is to check the T-spline rules and if it's needed
            %we must creat a new edge
            %The default code in this section is premitive it can be
            %improved by findobj methods
            index=edgeOrientparameter.Number;
            switch direct
                case 'horizontal'
                    Obstacle=0;
                    for j=index-1:-1:1
                        group=findobj(obj.edge,'Direction','horizontal','DirectionParameterCoordinate',obj.U(j));
                        if isempty(group)
                            continue
                        else
                            for c=1:max(size(group))
                                if group(c).FirstTip.Number<=newknot.Number && group(c).LastTip.Number>=newknot.Number
                                    Obstacle=1;
                                    if ~isempty(findobj(obj.controlpoint,'CPCV',newknot,'CPCU',obj.U(j)))
                                        
                                        sideControlPoint=findobj(obj.controlpoint,'CPCV',newknot,'CPCU',obj.U(j));
                                        sideEdge=Edge('vertical');
                                        
                                        sideEdge.FirstControlPoint=sideControlPoint;
                                        sideEdge.LastControlPoint=newControlPoint;
                                        sideEdge.FirstTip=obj.U(j);
                                        sideEdge.LastTip=edgeOrientparameter;
                                        sideEdge.DirectionParameterCoordinate=newknot;
                                        
                                        n=max(size(obj.edge));
                                        edgeArray=Edge.empty(n+1,0);
                                        edgeArray(1:n)=obj.edge;
                                        edgeArray(n+1)=sideEdge;
                                        obj.edge=edgeArray;
                                        
                                        sideControlPoint.DownPoint=newControlPoint;
                                        sideControlPoint.DownEdge=sideEdge;
                                        
                                        newControlPoint.UpPoint=sideControlPoint;
                                        newControlPoint.UpEdge=sideEdge;
                                    end
                                    break
                                end
                            end
                            if Obstacle==1
                                break
                            end
                        end
                    end
                    
                    Obstacle=0;
                    for j=index+1:max(size(obj.U))
                        group=findobj(obj.edge,'Direction','horizontal','DirectionParameterCoordinate',obj.U(j));
                        if isempty(group)
                            continue
                        else
                            for c=1:max(size(group))
                                if group(c).FirstTip.Number<=newknot.Number && group(c).LastTip.Number>=newknot.Number
                                    Obstacle=1;
                                    if ~isempty(findobj(obj.controlpoint,'CPCV',newknot,'CPCU',obj.U(j)))
                                        
                                        sideControlPoint=findobj(obj.controlpoint,'CPCV',newknot,'CPCU',obj.U(j));
                                        sideEdge=Edge('vertical');
                                        
                                        sideEdge.FirstControlPoint=newControlPoint;
                                        sideEdge.LastControlPoint=sideControlPoint;
                                        sideEdge.LastTip=obj.U(j);
                                        sideEdge.FirstTip=edgeOrientparameter;
                                        sideEdge.DirectionParameterCoordinate=newknot;
                                        
                                        n=max(size(obj.edge));
                                        edgeArray=Edge.empty(n+1,0);
                                        edgeArray(1:n)=obj.edge;
                                        edgeArray(n+1)=sideEdge;
                                        obj.edge=edgeArray;
                                        
                                        sideControlPoint.UpPoint=newControlPoint;
                                        sideControlPoint.UpEdge=sideEdge;
                                        
                                        newControlPoint.DownPoint=sideControlPoint;
                                        newControlPoint.DownEdge=sideEdge;
                                    end
                                    break
                                end
                            end
                            if Obstacle==1
                                break
                            end
                        end
                    end
                    
                    
                    
                case 'vertical'
                    Obstacle=0;
                    for j=index-1:-1:1
                        group=findobj(obj.edge,'Direction','vertical','DirectionParameterCoordinate',obj.V(j));
                        if isempty(group)
                            continue
                        else
                            for c=1:max(size(group))
                                if group(c).FirstTip.Number<=newknot.Number && group(c).LastTip.Number>=newknot.Number
                                    Obstacle=1;
                                    if ~isempty(findobj(obj.controlpoint,'CPCU',newknot,'CPCV',obj.V(j)))
                                        
                                        sideControlPoint=findobj(obj.controlpoint,'CPCU',newknot,'CPCV',obj.V(j));
                                        sideEdge=Edge('horizontal');
                                        
                                        sideEdge.FirstControlPoint=sideControlPoint;
                                        sideEdge.LastControlPoint=newControlPoint;
                                        sideEdge.FirstTip=obj.V(j);
                                        sideEdge.LastTip=edgeOrientparameter;
                                        sideEdge.DirectionParameterCoordinate=newknot;
                                        
                                        n=max(size(obj.edge));
                                        edgeArray=Edge.empty(n+1,0);
                                        edgeArray(1:n)=obj.edge;
                                        edgeArray(n+1)=sideEdge;
                                        obj.edge=edgeArray;
                                        
                                        sideControlPoint.RightPoint=newControlPoint;
                                        sideControlPoint.RightEdge=sideEdge;
                                        
                                        newControlPoint.LeftPoint=sideControlPoint;
                                        newControlPoint.LeftEdge=sideEdge;
                                    end
                                    break
                                end
                            end
                            if Obstacle==1
                                break
                            end
                        end
                    end
                    Obstacle=0;
                    for j=index+1:max(size(obj.V))
                        group=findobj(obj.edge,'Direction','vertical','DirectionParameterCoordinate',obj.V(j));
                        if isempty(group)
                            continue
                        else
                            for c=1:max(size(group))
                                if group(c).FirstTip.Number<=newknot.Number && group(c).LastTip.Number>=newknot.Number
                                    Obstacle=1;
                                    if ~isempty(findobj(obj.controlpoint,'CPCU',newknot,'CPCV',obj.V(j)))
                                        
                                        sideControlPoint=findobj(obj.controlpoint,'CPCU',newknot,'CPCV',obj.V(j));
                                        sideEdge=Edge('horizontal');
                                        
                                        sideEdge.FirstControlPoint=newControlPoint;
                                        sideEdge.LastControlPoint=sideControlPoint;
                                        sideEdge.LastTip=obj.V(j);
                                        sideEdge.FirstTip=edgeOrientparameter;
                                        sideEdge.DirectionParameterCoordinate=newknot;
                                        
                                        n=max(size(obj.edge));
                                        edgeArray=Edge.empty(n+1,0);
                                        edgeArray(1:n)=obj.edge;
                                        edgeArray(n+1)=sideEdge;
                                        obj.edge=edgeArray;
                                        
                                        
                                        sideControlPoint.LeftPoint=newControlPoint;
                                        sideControlPoint.LeftEdge=sideEdge;
                                        
                                        newControlPoint.RightPoint=sideControlPoint;
                                        newControlPoint.RightEdge=sideEdge;
                                    end
                                    break
                                end
                            end
                            if Obstacle==1
                                break
                            end
                        end
                    end
            end
            for i=1:numel(obj.controlpoint)
                obj.controlpoint(i).Number=i;
            end
            for i=1:numel(obj.edge)
                obj.edge(i).Number=i;
            end
            
        end
        %%
        function dataPoints = extract_Data(obj)
            resolution = 100;
            step=1/(resolution-1);
            uu = linspace(obj.U(1).Value,obj.U(end).Value-step,resolution);
            vv = linspace(obj.V(1).Value,obj.V(end).Value-step,resolution);
            dataPoints = zeros(numel(uu),numel(vv),3);
            
            for i = 1:numel(uu)
                for j = 1:numel(vv)
                    try
                        [x,y,z] = Calculate(obj,uu(i),vv(j));
                        dataPoints(i,j,1)=x; dataPoints(i,j,2)=y; dataPoints(i,j,3)=z;
                    catch
                        error(['i = ',num2str(i),', j = ',num2str(j)])
                    end
                end
            end
        end
        %%
        function obj=LocalInsertPoint(obj,edgenumber,varargin)% this function is used to add control points by the local refinement algorithm which is introduced in 2004 by sederberg
            %%
            %there are three types of using this function
            %% FIRST %%
            % obj = InsertPoint(obj,edgenumber)
            % to simplified the using of inserting control point user can
            % add T-spline object and the number of edge that user wants to
            % split then add a control point. parameter of new control point
            % is in the middle of former(selected) edge
            %% SECOND %%
            % obj = InsertPoint(obj,edgenumber,knot)
            % the inputs are T-spline object, the number of edge that user
            % wants to split and a existed knot which it's parameter is
            % locatad on edge.
            % this type of inputs is mostly used by program and functions in some
            % situations
            %% THIRD %%
            % obj = InsertPoint(obj,edgenumber,value)
            % the inputs are like previous type except value which is the
            % parameter value
            %this type mostly used by other programs and functions
            %%
            edgepopulation=numel(obj.edge);
            %to find the edge(the variable 'i' is important that indicate
            %the index of edge number
            for i=1:edgepopulation
                if obj.edge(i).Number==edgenumber
                    SelectedEdge=obj.edge(i);
                    break
                end
            end
            %
            edgeOrientparameter=SelectedEdge.DirectionParameterCoordinate;
            direct=SelectedEdge.Direction;
            %this section will define domain of point that would be
            %inserted
            if nargin==2
                switch direct
                    case 'horizontal'
                        
                        Val = (SelectedEdge.FirstTip.Value + SelectedEdge.LastTip.Value )/2;
                        status = 'notFound';
                        for i =  SelectedEdge.FirstTip.Number+1 : SelectedEdge.LastTip.Number-1
                            
                            if (Val-obj.V(i).Value) < 0.000001
                                newknot = obj.V(i);
                                status = 'Found';
                                break
                            end
                            
                        end
                        if strcmp(status,'notFound')
                            newknot=Vknot;
                            insertAfter(newknot,SelectedEdge.FirstTip);
                            insertBefore(newknot,SelectedEdge.LastTip);
                            %resorting in V
                            c=1;
                            knot=obj.V(1);
                            vv=Vknot.empty(numel(obj.V)+1,0);
                            vv(1)=knot;
                            while(true)
                                knot=knot.NextKnot;
                                if c<=numel(obj.V)
                                    vv(c+1)=knot;
                                    knot.Number=c+1;
                                    c=c+1;
                                else
                                    break
                                end
                            end
                            obj.V=vv;
                            newknot.Value=Val;
                        end
                    case 'vertical'
                      Val = (SelectedEdge.FirstTip.Value + SelectedEdge.LastTip.Value )/2;
                        status = 'notFound';
                        for i =  SelectedEdge.FirstTip.Number+1 : SelectedEdge.LastTip.Number-1
                            
                            if (Val-obj.U(i).Value) < 0.000001
                                newknot = obj.U(i);
                                status = 'Found';
                                break
                            end
                            
                        end
                        if strcmp(status,'notFound')
                            newknot=Uknot;
                            insertAfter(newknot,SelectedEdge.FirstTip);
                            insertBefore(newknot,SelectedEdge.LastTip);
                            %resorting in V
                            c=1;
                            knot=obj.U(1);
                            uu=Uknot.empty(numel(obj.U)+1,0);
                            uu(1)=knot;
                            while(true)
                                knot=knot.NextKnot;
                                if c<=numel(obj.U)
                                    uu(c+1)=knot;
                                    knot.Number=c+1;
                                    c=c+1;
                                else
                                    break
                                end
                            end
                            obj.U=uu;
                            newknot.Value=Val;
                        end
                end
                
            elseif nargin==3
                newknot=varargin{1};%%IT MIGHT BE WRONG
            elseif nargin==4
                switch direct
                    case 'horizontal'
                        newknot=Vknot;
                        k=SelectedEdge.FirstTip.Number;
                        j=SelectedEdge.LastTip.Number;
                        value=varargin{2};
                        while(true)
                            if value<obj.V(k+1).Value && value>=obj.V(k).Value
                                newknot.Number=k+1;
                                newknot.Value=value;
                                v=Vknot.empty(numel(obj.V)+1,0);
                                v(1:k)=obj.V(1:k);
                                v(k+1)=newknot;
                                v(k+2:numel(obj.V)+1)=obj.V(k+1:end);
                                obj.V=v;
                                break
                            end
                            k=k+1;
                        end
                        insertAfter(newknot,obj.V(k));
                        insertBefore(newknot,obj.V(k+1));
                        %resorting in V
                        for c=1:numel(obj.V)
                            obj.V(c).Number=c;
                        end
                        
                    case 'vertical'
                        newknot=Uknot;
                        k=SelectedEdge.FirstTip.Number;
                        j=SelectedEdge.LastTip.Number;
                        value=varargin{2};
                        while(true)
                            if value<obj.U(k+1).Value && value>=obj.U(k).Value
                                newknot.Number=k+1;
                                newknot.Value=value;
                                u=Uknot.empty(numel(obj.U)+1,0);
                                u(1:k)=obj.U(1:k);
                                u(k+1)=newknot;
                                u(k+2:numel(obj.U)+1)=obj.U(k+1:end);
                                obj.U=u;
                                break
                            end
                            k=k+1;
                        end
                        insertAfter(newknot,obj.U(k));
                        insertBefore(newknot,obj.U(k+1));
                        %resorting in U
                        for c=1:numel(obj.U)
                            obj.U(c).Number=c;
                        end
                end
            end
            %%
            
            %%
            if strcmp(direct,'horizontal')
                newCP=ControlPoint(0,0,0,0,SelectedEdge.FirstControlPoint.Uvector,[SelectedEdge.FirstControlPoint.Vvector(2) SelectedEdge.FirstControlPoint.CPCV newknot SelectedEdge.LastControlPoint.CPCV SelectedEdge.LastControlPoint.Vvector(4)]);
                newCP.Number=numel(obj.controlpoint)+1;
                obj.controlpoint=[obj.controlpoint,newCP];
                oldCP=ControlPoint.empty; newEQ={};
                %
                if  ~isempty(SelectedEdge.FirstControlPoint.LeftPoint)
                    oldCP(end+1)=SelectedEdge.FirstControlPoint.LeftPoint;
                    [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'horizontal');
                    newEQ{end+1,1}=A; newEQ{end,2}=B;
                end
                %
                oldCP(end+1)=SelectedEdge.FirstControlPoint;
                [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'horizontal');
                newEQ{end+1,1}=A; newEQ{end,2}=B;
                %
                oldCP(end+1)=SelectedEdge.LastControlPoint;
                [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'horizontal');
                newEQ{end+1,1}=A; newEQ{end,2}=B;
                %
                if  ~isempty(SelectedEdge.LastControlPoint.RightPoint)
                    oldCP(end+1)=SelectedEdge.LastControlPoint.RightPoint;
                    [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'horizontal');
                    newEQ{end+1,1}=A; newEQ{end,2}=B;
                end
                %%
                newEdge=Edge('horizontal'); newEdge.LastTip=SelectedEdge.LastTip; newEdge.FirstTip=newknot; newEdge.LastControlPoint=SelectedEdge.LastControlPoint;newEdge.FirstControlPoint=newCP;
                newEdge.DirectionParameterCoordinate=SelectedEdge.DirectionParameterCoordinate;
                
                SelectedEdge.LastTip=newknot; SelectedEdge.LastControlPoint=newCP; newEdge.FirstControlPoint=newCP;
                SelectedEdge.FirstControlPoint.RightPoint=newCP; SelectedEdge.LastControlPoint=newCP;
                obj.edge=[obj.edge,newEdge];
                newEdge.Number=numel(obj.edge);
                %%
                newCP.LeftEdge=SelectedEdge;
                newCP.LeftPoint=SelectedEdge.FirstControlPoint;
                newCP.RightPoint=newEdge.LastControlPoint;
                newCP.RightEdge=newEdge;
                newCP.Uvector=newCP.LeftPoint.Uvector; %fake vector for newCP
                newCP.Vvector=[newCP.LeftPoint.Vvector(2),newCP.LeftPoint.CPCV,newknot,newCP.RightPoint.CPCV,newCP.RightPoint.Vvector(4)]; %fake vector
                newCP.LeftEdge=SelectedEdge;
                %%
                newEdge.LastControlPoint.LeftPoint=newCP;
                newEdge.LastControlPoint.LeftEdge=newEdge;
                %%
                obj=Remesh(obj,newCP);
                %%
                
            elseif strcmp(direct,'vertical')
                newCP=ControlPoint(0,0,0,0,[SelectedEdge.FirstControlPoint.Uvector(2) SelectedEdge.FirstControlPoint.CPCU newknot SelectedEdge.LastControlPoint.CPCU SelectedEdge.LastControlPoint.Uvector(4)],SelectedEdge.FirstControlPoint.Vvector);
                newCP.Number=numel(obj.controlpoint)+1;
                obj.controlpoint= [obj.controlpoint,newCP];
                oldCP=ControlPoint.empty; newEQ={};
                %
                if  ~isempty(SelectedEdge.FirstControlPoint.UpPoint)
                    oldCP(end+1)=SelectedEdge.FirstControlPoint.UpPoint;
                    [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'vertical');
                    newEQ{end+1,1}=A; newEQ{end,2}=B;
                end
                %
                oldCP(end+1)=SelectedEdge.FirstControlPoint;
                [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'vertical');
                newEQ{end+1,1}=A; newEQ{end,2}=B;
                %
                oldCP(end+1)=SelectedEdge.LastControlPoint;
                [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'vertical');
                newEQ{end+1,1}=A; newEQ{end,2}=B;
                %
                if  ~isempty(SelectedEdge.LastControlPoint.DownPoint)
                    oldCP(end+1)=SelectedEdge.LastControlPoint.DownPoint;
                    [A,B]=Coef(obj,oldCP(end),oldCP(end).Uvector,oldCP(end).Vvector,newknot,'vertical');
                    newEQ{end+1,1}=A; newEQ{end,2}=B;
                end
                %%
                newEdge=Edge('vertical'); newEdge.LastTip=SelectedEdge.LastTip; newEdge.FirstTip=newknot; newEdge.LastControlPoint=SelectedEdge.LastControlPoint;newEdge.FirstControlPoint=newCP;
                newEdge.DirectionParameterCoordinate=SelectedEdge.DirectionParameterCoordinate;
                SelectedEdge.LastTip=newknot; SelectedEdge.LastControlPoint=newCP; newEdge.FirstControlPoint=newCP;
                SelectedEdge.FirstControlPoint.DownPoint=newCP; SelectedEdge.LastControlPoint=newCP;
                obj.edge=[obj.edge,newEdge];
                newEdge.Number=numel(obj.edge);
                %%
                newCP.UpEdge=SelectedEdge; newCP.UpPoint=SelectedEdge.FirstControlPoint; newCP.DownEdge=newEdge; newCP.DownPoint=newEdge.LastControlPoint;
                
                %%
                newEdge.LastControlPoint.UpPoint=newCP;
                newEdge.LastControlPoint.UpEdge=newEdge;
                %%
                newCP.DownPoint=newEdge.LastControlPoint;
                newCP.UpPoint=SelectedEdge.FirstControlPoint;
                newCP.DownEdge=newEdge;
                newCP.Vvector=newCP.UpPoint.Vvector; %fake vector for newCP
                newCP.Uvector=[newCP.UpPoint.Uvector(2),newCP.UpPoint.CPCU,newknot,newCP.DownPoint.CPCU,newCP.DownPoint.Uvector(4)]; %fake vector
                newCP.UpEdge=SelectedEdge;
                %%
                obj=Remesh(obj,newCP);
                %%
            end
            %%
            newPoints=ControlPoint.empty;
            newPoints(end+1) = newCP;
            while(true)
                status='finished';
                for i=1:numel(oldCP)
                    for j=1:size(newEQ,2)
                        if ~isempty(newEQ{i,j})
                            %%%%%%%%%%%%%% vertical check for each control
                            %%%%%%%%%%%%%% point
                            for L=2:-1:1
                                [umesh,vmesh]=CheckDomain(obj,newEQ{i,j}{4}.Number);
                                uu=newEQ{i,j}{1}; vv=newEQ{i,j}{2};
                                %%
                                if umesh(L).Number<uu(L).Number
                                    Point=ControlPoint(0,0,0,0,[umesh(1),umesh(2),uu(L),umesh(3),umesh(4)],newEQ{i,j}{4}.Vvector);
                                    Point.Number=numel(obj.controlpoint)+1;
                                    obj.controlpoint= [obj.controlpoint,Point];
                                    Point.Vvector=newEQ{i,j}{4}.Vvector;
                                    newPoints(end+1)=Point;
                                    %%
                                    if L==2
                                        Point.Uvector=[umesh(1),umesh(2),uu(L),umesh(3),umesh(4)];
                                        obj=Remesh(obj,Point);
                                    else
                                        Point.Uvector=[umesh(1),umesh(1),uu(L),umesh(2),umesh(3)]; %fake vector
                                        obj=Remesh(obj,Point);
                                        [u,v]=CheckDomain(obj,Point.Number);
                                        Point.Uvector=u;
                                        Point.Vvector=v;
                                    end
                                    if ~isempty(Point.LeftPoint) && isempty(find(Point.LeftPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.LeftPoint==newPoints))
                                        [A,B]=Coef(obj,Point.LeftPoint,Point.LeftPoint.Uvector,Point.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.LeftPoint==oldCP)) 
                                            oldCP(end+1)=Point.LeftPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.LeftPoint.LeftPoint) && isempty(find(Point.LeftPoint.LeftPoint==newPoints))
                                            [A,B]=Coef(obj,Point.LeftPoint.LeftPoint,Point.LeftPoint.LeftPoint.Uvector,Point.LeftPoint.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.LeftPoint.LeftPoint==oldCP))
                                                oldCP(end+1)=Point.LeftPoint.LeftPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.RightPoint) && isempty(find(Point.RightPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.RightPoint==newPoints))
                                        [A,B]=Coef(obj,Point.RightPoint,Point.RightPoint.Uvector,Point.RightPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.RightPoint==oldCP))
                                            oldCP(end+1)=Point.RightPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.RightPoint.RightPoint) && isempty(find(Point.RightPoint.RightPoint==newPoints))
                                            [A,B]=Coef(obj,Point.RightPoint.RightPoint,Point.RightPoint.RightPoint.Uvector,Point.RightPoint.RightPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.RightPoint.RightPoint==oldCP))
                                                oldCP(end+1)=Point.RightPoint.RightPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%new update
                                    if ~isempty(Point.UpPoint) && isempty(find(Point.UpPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.UpPoint==newPoints))
                                        [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.Uvector,Point.UpPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.UpPoint==oldCP))
                                            oldCP(end+1)=Point.UpPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.UpPoint.UpPoint) && isempty(find(Point.UpPoint.UpPoint==newPoints))
                                            [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.UpPoint.Uvector,Point.UpPoint.UpPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.UpPoint.UpPoint==oldCP))
                                                oldCP(end+1)=Point.UpPoint.UpPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.DownPoint) && isempty(find(Point.DownPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.DownPoint==newPoints))
                                        [A,B]=Coef(obj,Point.DownPoint,Point.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.DownPoint==oldCP))
                                            oldCP(end+1)=Point.DownPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.DownPoint.DownPoint) && isempty(find(Point.DownPoint.DownPoint==newPoints))
                                            [A,B]=Coef(obj,Point.DownPoint.DownPoint,Point.DownPoint.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.DownPoint.DownPoint==oldCP))
                                                oldCP(end+1)=Point.DownPoint.DownPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%end update
                                    
                                    status='unfinished';
                                elseif umesh(L).Number>uu(L).Number
                                    [A,B]=Coef(obj,newEQ{i,j}{4},uu,newEQ{i,j}{2},umesh(L),'vertical');
                                    A{3}=A{3}*newEQ{i,j}{3}; B{3}=B{3}*newEQ{i,j}{3};
                                    newEQ{i,j}=B; c=0;
                                    for k=j+1:size(newEQ,2)
                                        if isempty(newEQ{i,k})
                                            newEQ{i,k}=A; c=1;
                                            break
                                        end
                                    end
                                    if c==0
                                        newEQ{i,end+1}=A;
                                    end
                                    status='unfinished';
                                end
                                
                            end
                            for L=4:1:5
                                [umesh,vmesh]=CheckDomain(obj,newEQ{i,j}{4}.Number);
                                uu=newEQ{i,j}{1}; vv=newEQ{i,j}{2};
                                %%
                                if umesh(L).Number>uu(L).Number
                                    Point=ControlPoint(0,0,0,0,[umesh(2),umesh(3),uu(L),umesh(4),umesh(5)],newEQ{i,j}{4}.Vvector);
                                    Point.Number=numel(obj.controlpoint)+1;
                                    obj.controlpoint=[obj.controlpoint,Point];
                                    Point.Vvector=newEQ{i,j}{4}.Vvector;
                                    newPoints(end+1)=Point;
                                    if L==4
                                        Point.Uvector=[umesh(2),umesh(3),uu(L),umesh(4),umesh(5)];
                                        obj=Remesh(obj,Point);
                                    else
                                        Point.Uvector=[umesh(3),umesh(4),uu(L),umesh(5),umesh(5)];%fake vector
                                        obj=Remesh(obj,Point);
                                        [u,v]=CheckDomain(obj,Point.Number);
                                        Point.Uvector=u;
                                        Point.Vvector=v;
                                        %%
                                    end
                                    if ~isempty(Point.LeftPoint) && isempty(find(Point.LeftPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.LeftPoint==newPoints))
                                        [A,B]=Coef(obj,Point.LeftPoint,Point.LeftPoint.Uvector,Point.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.LeftPoint==oldCP))
                                            oldCP(end+1)=Point.LeftPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.LeftPoint.LeftPoint) && isempty(find(Point.LeftPoint.LeftPoint==newPoints))
                                            [A,B]=Coef(obj,Point.LeftPoint.LeftPoint,Point.LeftPoint.LeftPoint.Uvector,Point.LeftPoint.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.LeftPoint.LeftPoint==oldCP))
                                                oldCP(end+1)=Point.LeftPoint.LeftPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.RightPoint) && isempty(find(Point.RightPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.RightPoint==newPoints))
                                        [A,B]=Coef(obj,Point.RightPoint,Point.RightPoint.Uvector,Point.RightPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.RightPoint==oldCP))
                                            oldCP(end+1)=Point.RightPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.RightPoint.RightPoint) && isempty(find(Point.RightPoint.RightPoint==newPoints))
                                            [A,B]=Coef(obj,Point.RightPoint.RightPoint,Point.RightPoint.RightPoint.Uvector,Point.RightPoint.RightPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.RightPoint.RightPoint==oldCP))
                                                oldCP(end+1)=Point.RightPoint.RightPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%new update
                                    if ~isempty(Point.UpPoint) && isempty(find(Point.UpPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.UpPoint==newPoints))
                                        [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.Uvector,Point.UpPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.UpPoint==oldCP))
                                            oldCP(end+1)=Point.UpPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.UpPoint.UpPoint) && isempty(find(Point.UpPoint.UpPoint==newPoints))
                                            [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.UpPoint.Uvector,Point.UpPoint.UpPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.UpPoint.UpPoint==oldCP))
                                                oldCP(end+1)=Point.UpPoint.UpPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.DownPoint) && isempty(find(Point.DownPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.DownPoint==newPoints))
                                        [A,B]=Coef(obj,Point.DownPoint,Point.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.DownPoint==oldCP))
                                            oldCP(end+1)=Point.DownPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.DownPoint.DownPoint) && isempty(find(Point.DownPoint.DownPoint==newPoints,1))
                                            [A,B]=Coef(obj,Point.DownPoint.DownPoint,Point.DownPoint.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.DownPoint.DownPoint==oldCP,1))
                                                oldCP(end+1)=Point.DownPoint.DownPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%end update
                                    status='unfinished';
                                elseif umesh(L).Number<uu(L).Number
                                    [A,B]=Coef(obj,newEQ{i,j}{4},uu,newEQ{i,j}{2},umesh(L),'vertical');
                                    A{3}=A{3}*newEQ{i,j}{3}; B{3}=B{3}*newEQ{i,j}{3};
                                    newEQ{i,j}=A; c=0;
                                    for k=j+1:size(newEQ,2)
                                        if isempty(newEQ{i,k})
                                            newEQ{i,k}=B; c=1;
                                            break
                                        end
                                    end
                                    if c==0
                                        newEQ{i,end+1}=B;
                                    end
                                    status='unfinished';
                                end
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%% horizontal check for
                            %%%%%%%%%%%%%%%%%%%%%%%%% each control point
                            for L=2:-1:1
                                [umesh,vmesh]=CheckDomain(obj,newEQ{i,j}{4}.Number);
                                uu=newEQ{i,j}{1}; vv=newEQ{i,j}{2};
                                %%
                                if vmesh(L).Number<vv(L).Number
                                    Point=ControlPoint(0,0,0,0,newEQ{i,j}{4}.Uvector,[vmesh(1),vmesh(2),vv(L),vmesh(3),vmesh(4)]);
                                    Point.Number=numel(obj.controlpoint)+1;
                                    obj.controlpoint= [obj.controlpoint,Point];
                                    Point.Vvector=newEQ{i,j}{4}.Vvector; %% edges is not remembered
                                    newPoints(end+1)=Point;
                                    if L==2
                                        Point.Vvector=[vmesh(1),vmesh(2),vv(L),vmesh(3),vmesh(4)];
                                        obj=Remesh(obj,Point);
                                    else
                                        Point.Vvector=[vmesh(1),vmesh(1),vv(L),vmesh(2),vmesh(3)];
                                        obj=Remesh(obj,Point);
                                        [u,v]=CheckDomain(obj,Point.Number);
                                        Point.Vvector=v;
                                        Point.Uvector=u;
                                    end
                                    %%
                                    if ~isempty(Point.UpPoint) && isempty(find(Point.UpPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.UpPoint==newPoints))
                                        [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.Uvector,Point.UpPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.UpPoint==oldCP))
                                            oldCP(end+1)=Point.UpPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.UpPoint.UpPoint) && isempty(find(Point.UpPoint.UpPoint==newPoints))
                                            [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.UpPoint.Uvector,Point.UpPoint.UpPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.UpPoint.UpPoint==oldCP))
                                                oldCP(end+1)=Point.UpPoint.UpPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.DownPoint) && isempty(find(Point.DownPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.DownPoint==newPoints))
                                        [A,B]=Coef(obj,Point.DownPoint,Point.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.DownPoint==oldCP))
                                            oldCP(end+1)=Point.DownPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.DownPoint.DownPoint) && isempty(find(Point.DownPoint.DownPoint==newPoints))
                                            [A,B]=Coef(obj,Point.DownPoint.DownPoint,Point.DownPoint.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.DownPoint.DownPoint==oldCP))
                                                oldCP(end+1)=Point.DownPoint.DownPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%new update
                                    if ~isempty(Point.LeftPoint) && isempty(find(Point.LeftPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.LeftPoint==newPoints))
                                        [A,B]=Coef(obj,Point.LeftPoint,Point.LeftPoint.Uvector,Point.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.LeftPoint==oldCP))
                                            oldCP(end+1)=Point.LeftPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.LeftPoint.LeftPoint) && isempty(find(Point.LeftPoint.LeftPoint==newPoints))
                                            [A,B]=Coef(obj,Point.LeftPoint.LeftPoint,Point.LeftPoint.LeftPoint.Uvector,Point.LeftPoint.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.LeftPoint.LeftPoint==oldCP))
                                                oldCP(end+1)=Point.LeftPoint.LeftPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.RightPoint) && isempty(find(Point.RightPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.RightPoint==newPoints))
                                        [A,B]=Coef(obj,Point.RightPoint,Point.RightPoint.Uvector,Point.RightPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.RightPoint==oldCP))
                                            oldCP(end+1)=Point.RightPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.RightPoint.RightPoint) && isempty(find(Point.RightPoint.RightPoint==newPoints))
                                            [A,B]=Coef(obj,Point.RightPoint.RightPoint,Point.RightPoint.RightPoint.Uvector,Point.RightPoint.RightPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.RightPoint.RightPoint==oldCP))
                                                oldCP(end+1)=Point.RightPoint.RightPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%end update
                                    status='unfinished';
                                    %%
                                elseif vmesh(L).Number>vv(L).Number
                                    [A,B]=Coef(obj,newEQ{i,j}{4},newEQ{i,j}{1},vv,vmesh(L),'horizontal');
                                    A{3}=A{3}*newEQ{i,j}{3}; B{3}=B{3}*newEQ{i,j}{3};
                                    newEQ{i,j}=B; c=0;
                                    for k=j+1:size(newEQ,2)
                                        if isempty(newEQ{i,k})
                                            newEQ{i,k}=A; c=1;
                                            break
                                        end
                                    end
                                    if c==0
                                        newEQ{i,end+1}=A;
                                    end
                                    status='unfinished';
                                end
                                
                            end
                            for L=4:1:5
                                [umesh,vmesh]=CheckDomain(obj,newEQ{i,j}{4}.Number);
                                uu=newEQ{i,j}{1}; vv=newEQ{i,j}{2};
                                %%
                                if vmesh(L).Number>vv(L).Number
                                    Point=ControlPoint(0,0,0,0,newEQ{i,j}{4}.Uvector,[vmesh(2),vmesh(3),vv(L),vmesh(4),vmesh(5)]);
                                    Point.Number=numel(obj.controlpoint)+1;
                                    obj.controlpoint= [obj.controlpoint,Point];
                                    Point.Uvector=newEQ{i,j}{4}.Uvector; %% edges is not remembered
                                    newPoints(end+1) = Point;
                                    if L==4
                                        Point.Vvector=[vmesh(2),vmesh(3),vv(L),vmesh(4),vmesh(5)];
                                        obj=Remesh(obj,Point);
                                    else
                                        Point.Vvector=[vmesh(3),vmesh(4),vv(L),vmesh(5),vmesh(5)];
                                        obj=Remesh(obj,Point);
                                        [u,v]=CheckDomain(obj,Point.Number);
                                        Point.Vvector=v;
                                        Point.Uvector=u;
                                    end
                                    %%
                                    if ~isempty(Point.UpPoint) && isempty(find(Point.UpPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.UpPoint==newPoints))
                                        [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.Uvector,Point.UpPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.UpPoint==oldCP))
                                            oldCP(end+1)=Point.UpPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.UpPoint.UpPoint) && isempty(find(Point.UpPoint.UpPoint==newPoints))
                                            [A,B]=Coef(obj,Point.UpPoint,Point.UpPoint.UpPoint.Uvector,Point.UpPoint.UpPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.UpPoint.UpPoint==oldCP))
                                                oldCP(end+1)=Point.UpPoint.UpPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.DownPoint) && isempty(find(Point.DownPoint.Uvector == Point.CPCU, 1)) && isempty(find(Point.DownPoint==newPoints))
                                        [A,B]=Coef(obj,Point.DownPoint,Point.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                        if isempty(find(Point.DownPoint==oldCP))
                                            oldCP(end+1)=Point.DownPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.DownPoint.DownPoint) && isempty(find(Point.DownPoint.DownPoint==newPoints))
                                            [A,B]=Coef(obj,Point.DownPoint.DownPoint,Point.DownPoint.DownPoint.Uvector,Point.DownPoint.Vvector,Point.CPCU,'vertical');
                                            if isempty(find(Point.DownPoint.DownPoint==oldCP))
                                                oldCP(end+1)=Point.DownPoint.DownPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%new update
                                    if ~isempty(Point.LeftPoint) && isempty(find(Point.LeftPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.LeftPoint==newPoints))
                                        [A,B]=Coef(obj,Point.LeftPoint,Point.LeftPoint.Uvector,Point.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.LeftPoint==oldCP))
                                            oldCP(end+1)=Point.LeftPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.LeftPoint.LeftPoint) && isempty(find(Point.LeftPoint.LeftPoint==newPoints))
                                            [A,B]=Coef(obj,Point.LeftPoint.LeftPoint,Point.LeftPoint.LeftPoint.Uvector,Point.LeftPoint.LeftPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.LeftPoint.LeftPoint==oldCP))
                                                oldCP(end+1)=Point.LeftPoint.LeftPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    if ~isempty(Point.RightPoint) && isempty(find(Point.RightPoint.Vvector == Point.CPCV, 1)) && isempty(find(Point.RightPoint==newPoints))
                                        [A,B]=Coef(obj,Point.RightPoint,Point.RightPoint.Uvector,Point.RightPoint.Vvector,Point.CPCV,'horizontal');
                                        if isempty(find(Point.RightPoint==oldCP))
                                            oldCP(end+1)=Point.RightPoint;
                                            newEQ{end+1,1}=A; newEQ{end,2}=B;
                                        end
                                        if ~isempty(Point.RightPoint.RightPoint) && isempty(find(Point.RightPoint.RightPoint==newPoints))
                                            [A,B]=Coef(obj,Point.RightPoint.RightPoint,Point.RightPoint.RightPoint.Uvector,Point.RightPoint.RightPoint.Vvector,Point.CPCV,'horizontal');
                                            if isempty(find(Point.RightPoint.RightPoint==oldCP))
                                                oldCP(end+1)=Point.RightPoint.RightPoint;
                                                newEQ{end+1,1}=A; newEQ{end,2}=B;
                                            end
                                        end
                                    end
                                    %%end update
                                    status='unfinished';
                                    %%
                                elseif vmesh(L).Number<vv(L).Number
                                    [A,B]=Coef(obj,newEQ{i,j}{4},newEQ{i,j}{1},vv,vmesh(L),'horizontal');
                                    A{3}=A{3}*newEQ{i,j}{3}; B{3}=B{3}*newEQ{i,j}{3};
                                    newEQ{i,j}=A; c=0;
                                    for k=j+1:size(newEQ,2)
                                        if isempty(newEQ{i,k})
                                            newEQ{i,k}=B; c=1;
                                            break
                                        end
                                    end
                                    if c==0
                                        newEQ{i,end+1}=B;
                                    end
                                    status='unfinished';
                                end
                            end
                        end
                    end
                end
                %%
                if strcmp(status,'finished')
                    break
                end
            end
            %% calculating control points
            for i=1:size(newEQ,1)
                for j=1:size(newEQ,2)
                    if ~isempty(newEQ{i,j})
                        newEQ{i,j}{end+1}=oldCP(i);
                    end
                end
            end
            
            CPs={}; CPs{end+1}=newEQ{1}{4};
            for i=1:numel(newEQ)
                if isempty(newEQ{i})
                    continue
                end
                c=0;
                for j=1:numel(CPs)
                    if CPs{j}==newEQ{i}{4}
                        c=1;
                    end
                end
                if c==0
                    CPs{end+1}=newEQ{i}{4};
                end
            end
            %
            coordinate=[];
            for i=1:numel(CPs)
                xw=0;    yw=0;    zw=0;    w=0;
                existance = 'No';
                for k = 1:numel(oldCP)
                    if CPs{i} == oldCP(k)
                        existance = 'Yes';
                        break
                    end
                end
                if strcmp(existance,'No')
                    xw = CPs{i}.XW;
                    yw = CPs{i}.YW;
                    zw = CPs{i}.ZW;
                    w = CPs{i}.W;
                end
                for j=1:numel(newEQ)
                    if isempty(newEQ{j})
                        continue
                    end
                    if CPs{i}==newEQ{j}{4}
                        xw=xw+newEQ{j}{5}.XW*newEQ{j}{3};
                        yw=yw+newEQ{j}{5}.YW*newEQ{j}{3};
                        zw=zw+newEQ{j}{5}.ZW*newEQ{j}{3};
                        w=w+newEQ{j}{5}.W*newEQ{j}{3};
                        newEQ{j}{4}.Uvector=newEQ{j}{1};
                        newEQ{j}{4}.Vvector=newEQ{j}{2};
                    end
                end
                coordinate(end+1,1)=xw;
                coordinate(end,2)=yw;
                coordinate(end,3)=zw;
                coordinate(end,4)=w;
                
            end
            
            for i=1:numel(CPs)
                CPs{i}.X=coordinate(i,1)/coordinate(i,4);
                CPs{i}.Y=coordinate(i,2)/coordinate(i,4);
                CPs{i}.Z=coordinate(i,3)/coordinate(i,4);
                CPs{i}.W=coordinate(i,4);
            end
            
        end
        %%
        function newobj=Merge(obj,number1,number2,normal1_number,obj2,number3,number4,normal2_number,continuity)%%%obj:first surface obj2:second surface
            %%%number1:first controlpoint located at begginning of the shared border in first surface; number2:second controlpoint located at end of the shared border in first surface
            %%%number3:first controlpoint located at begginning of the shared border in second surface; number4:second controlpoint located at end of the shared border in second surface
            %%%normal1: controlpoint located behind the shared border in first surface; normal2: controlpoint located behind of the shared border in second surface
            %%%continuity: C0 or C1 or C2
            
            %%to find controlpoints on surfaces
            point1=findobj(obj.controlpoint,'Number',number1);
            point2=findobj(obj.controlpoint,'Number',number2);
            normal1=findobj(obj.controlpoint,'Number',normal1_number);
            
            point3=findobj(obj2.controlpoint,'Number',number3);
            point4=findobj(obj2.controlpoint,'Number',number4);
            normal2=findobj(obj2.controlpoint,'Number',normal2_number);
            
            clear number1 number2 normal1_number number3 number4 normal2_number
            
            %%
            %in this section we try to make a order for surfaces it
            %means their vectors or grid become justify
            param1={point1.CPCU point1.CPCV};
            param2={point2.CPCU point2.CPCV};
            
            param3={point3.CPCU point3.CPCV};
            param4={point4.CPCU point4.CPCV};
            %%
            if param1{1}==param2{1}
                obj=SwitchVectors(obj);
            end
            Axial1='V';
            Transverse1='U';
            isoparm1=point1.CPCV;
            %%
            if param3{1}==param4{1}
                obj2=SwitchVectors(obj2);
            end
            Axial2='V';
            Transverse2='U';
            isoparm2=point3.CPCV;
            %%
            
            dnumber1=point2.CPCU.Number-point1.CPCU.Number;
            dnumber2=point4.CPCU.Number-point3.CPCU.Number;
            
            if dnumber1<0
                obj=ReverseKnotVector(obj,Transverse1);
            end
            if dnumber2<0
                obj2=ReverseKnotVector(obj2,Transverse2);
            end
            %%
            dnumNormal1=normal1.CPCV.Number-point1.CPCV.Number;
            if dnumNormal1>0
                obj=ReverseKnotVector(obj,Axial1);
            end
            
            dnumNormal2=normal2.CPCV.Number-point3.CPCV.Number;
            if dnumNormal2<0
                obj2=ReverseKnotVector(obj2,Axial2);
            end
            clear param1 param2 dnumber1 dnumber2 dnumNormal1 dnumNormal2
            %% finding knots of points in borders
            %for surface1
            sharevector1={};
            i=1;
            point=point1;
            while(true)
                P1{i}=point;
                sharevector1{i}=point.CPCU;
                if point==point2
                    break
                end
                i=i+1;
                point=point.DownPoint;
            end
            
            %for surface 2
            sharevector2={};
            i=1;
            point=point3;
            while(true)
                P2{i}=point;
                sharevector2{i}=point.CPCU;
                if point==point4
                    break
                end
                i=i+1;
                point=point.DownPoint;
            end
            %% finding intervals of borders and justify them for another surface
            Interval1=[];
            for i=1:length(sharevector1)-1
                Interval1(i)=sharevector1{i+1}.Value-sharevector1{i}.Value;
            end
            Interval2=[];
            for i=1:length(sharevector2)-1
                Interval2(i)=sharevector2{i+1}.Value-sharevector2{i}.Value;
            end
            %%
            ratio=sum(Interval1,2)/sum(Interval2,2);
            %%
            newInterval2=ratio*Interval2;
            %%
            localvalue1=[];
            localvalue1(1)=0;
            for i=1:length(Interval1)
                localvalue1(i+1)=localvalue1(i)+Interval1(i);
            end
            %%
            localvalue2=[];
            localvalue2(1)=0;
            for i=1:length(Interval2)
                localvalue2(i+1)=localvalue2(i)+newInterval2(i);
            end
            %%
            missed1in2=[]; missed2in1=[];
            %%a little bug we have here but we can ignore it
            L1=localvalue1(2:end-1);
            L2=localvalue2(2:end-1);
            
            status='notfound';
            i=1;c=0;
            while(true)
                for j=1:numel(L2)
                    if abs(L1(i)-L2(j))<0.00001
                        status='found';
                        if numel(L2)==1
                            c=1;
                            break
                        end
                        L2=L2(j+1:end);
                        break
                    end
                end
                if strcmp(status,'notfound')
                    missed1in2(end+1)=L1(i);
                end
                i=i+1;
                if i>numel(L1) || c==1
                    break
                end
                status='notfound';
            end
            
            L1=localvalue1(2:end-1);
            L2=localvalue2(2:end-1);
            status='notfound';
            i=1;c=0;
            while(true)
                for j=1:numel(L1)
                    if abs(L1(j)-L2(i))<0.00001
                        status='found';
                        if numel(L1)==1
                            c=1;
                            break
                        end
                        L1=L1(j+1:end);
                        break
                    end
                    
                end
                if strcmp(status,'notfound')
                    missed2in1(end+1)=L2(i);
                end
                i=i+1;
                if i>numel(L2) || c==1
                    break
                end
                status='notfound';
            end
            
            %%
            switch continuity
                case 'C0'
                    C=0;
                case 'C1'
                    C=1;
                case 'C2'
                    C=2;
            end
            %% inserting knots for surface 1
            point=point1;
            inknot1={};
            cp1={};
            Edge1={};
            for j=1:numel(missed2in1)
                
                value=point1.CPCU.Value+missed2in1(j);
                D=point.CPCV;
                
                groupEdge=findobj(obj.edge,'DirectionParameterCoordinate',D);
                for k=1:numel(groupEdge)
                    if abs(value-1)>0.00001
                        if groupEdge(k).FirstTip.Value<=value && groupEdge(k).LastTip.Value>value
                            Edge1{end+1}=[groupEdge(k).FirstTip,groupEdge(k).LastTip];
                            obj=InsertPoint(obj,groupEdge(k).Number,'Value',value);%%%%%%%%%%%%%%%%%%%%%%
                            inknot1{end+1}=obj.controlpoint(end).CPCU;
                            cp1{end+1}=obj.controlpoint(end);
                            break
                        end
                    else
                        if groupEdge(k).FirstTip.Value<=value && groupEdge(k).LastTip.Value>=value
                            Edge1{end+1}=[groupEdge(k).FirstTip,groupEdge(k).LastTip];
                            obj=InsertPoint(obj,groupEdge(k).Number,'Value',value);%%%%%%%%%%%%%%%%%%%%%%
                            inknot1{end+1}=obj.controlpoint(end).CPCU;
                            cp1{end+1}=obj.controlpoint(end);
                            break
                        end
                    end
                end
            end
            for i=1:numel(obj.controlpoint)
                obj.controlpoint(i).Number=i;
            end
            for i=1:numel(obj.edge)
                obj.edge(i).Number=i;
            end
            %% inserting knot in surface1 multiple times(C1 or C2)
            pp=point1;
            share1={};
            while(true)
                share1{end+1}=pp;
                if pp==point2
                    break
                else
                    pp=pp.DownPoint;
                end
            end
            
            for i=1:C
                knot=obj.V(share1{1}.CPCV.Number-i);
                obj=InsertRow(obj,knot);
                for j=1:numel(share1)
                    knotU=share1{j}.CPCU;
                    pp=findobj(obj.controlpoint,'CPCV',knot,'CPCU',knotU);
                    if isempty(pp)
                        ed=findobj(obj.edge,'DirectionParameterCoordinate',knot);
                        for k=1:numel(ed)
                            if knotU.Number>ed(k).FirstTip.Number && knotU.Number<ed(k).LastTip.Number
                                obj=InsertPoint(obj,ed(k).Number,knotU);
                                break
                            end
                        end
                    end
                end
            end
            
            %% inserting knots for surface 2
            point=point3;
            inknot2={};
            cp2={};
            Edge2={};
            for j=1:numel(missed1in2)
                
                value=point3.CPCU.Value+missed1in2(j)/ratio;
                D=findobj(obj2.V,'Number',point.CPCV.Number);
                
                groupEdge=findobj(obj2.edge,'DirectionParameterCoordinate',D);
                for k=1:numel(groupEdge)
                    if abs(value-1)>0.00001
                        if groupEdge(k).FirstTip.Value<=value && groupEdge(k).LastTip.Value>value
                            Edge2{end+1}=[groupEdge(k).FirstTip,groupEdge(k).LastTip];
                            obj2=InsertPoint(obj2,groupEdge(k).Number,'Value',value);%%%%%%%%%%%%%%%%%%%%%%
                            inknot2{end+1}=obj2.controlpoint(end).CPCU;
                            cp2{end+1}=obj2.controlpoint(end);
                            break
                        end
                    else
                        if groupEdge(k).FirstTip.Value<=value && groupEdge(k).LastTip.Value>=value
                            Edge2{end+1}=[groupEdge(k).FirstTip,groupEdge(k).LastTip];
                            obj2=InsertPoint(obj2,groupEdge(k).Number,'Value',value);%%%%%%%%%%%%%%%%%%%%%%
                            inknot2{end+1}=obj2.controlpoint(end).CPCU;
                            cp2{end+1}=obj2.controlpoint(end);
                            break
                        end
                    end
                end
            end
            for i=1:numel(obj2.controlpoint)
                obj2.controlpoint(i).Number=i;
            end
            for i=1:numel(obj2.edge)
                obj2.edge(i).Number=i;
            end
            %% inserting knot in surface2 multiple times(C1 or C2)%this has bug!!!!!!!!!!!!!!!!!!!
            
            pp=point3;
            share2={};
            while(true)
                share2{end+1}=pp;
                if pp==point4
                    break
                else
                    pp=pp.DownPoint;
                end
            end
            
            for i=1:C
                knot=obj2.V(share2{1}.CPCV.Number+i);
                obj2=InsertRow(obj2,knot);
                for j=1:numel(share2)
                    knotU=share2{j}.CPCU;
                    pp=findobj(obj.controlpoint,'CPCV',knot,'CPCU',knotU);
                    if isempty(pp)
                        ed=findobj(obj2.edge,'DirectionParameterCoordinate',knot);
                        for k=1:numel(ed)
                            if knotU.Number>ed(k).FirstTip.Number && knotU.Number<ed(k).LastTip.Number
                                obj2=InsertPoint(obj2,ed(k).Number,knotU);
                                break
                            end
                        end
                    end
                end
            end
            %% mixing two v knotvectors
            clear status param3 param4 k j c
            I1=obj.V;
            I2=obj2.V;
            interval=[];
            for i=1:numel(I2)-1
                interval(i)=I2(i+1).Value-I2(i).Value;
            end
            
            index1=isoparm1.Number;
            index2=isoparm2.Number;
            
            isoparm2.Value=isoparm1.Value;
            for i=index2-1:-1:1
                I2(i).Value=I2(i+1).Value-interval(i);
            end
            for i=index2+1:numel(I2)
                I2(i).Value=I2(i-1).Value+interval(i-1);
            end
            
            %%
            wholeknotV=Vknot.empty(numel(I1)+numel(I2),0);
            wholeknotV(1:numel(I1))=I1(1:end);
            wholeknotV(numel(I1)+1:numel(I1)+numel(I2))=I2(1:end);
            %% sorting of wholeknot
            i=2;
            while(true)
                if i>1 && wholeknotV(i-1).Value>wholeknotV(i).Value
                    pre=wholeknotV(i-1);
                    wholeknotV(i-1)=wholeknotV(i);
                    wholeknotV(i)=pre;
                    i=i-2;
                end
                i=i+1;
                if numel(wholeknotV)+1==i
                    break
                end
            end
            for i=1:numel(wholeknotV)
                wholeknotV(i).Number=i;
            end
            %% it was end of mixing two v vectors
            %% mixing two u vectors
            I1=obj.U; I2=obj2.U;
            knot1=point1.CPCU;  knot3=point3.CPCU;
            knot2=point2.CPCU;  knot4=point4.CPCU;
            
            %%
            for i=1:numel(I2)
                I2(i).Value=ratio*I2(i).Value;
            end
            %%
            index1=knot1.Number;
            index2=knot2.Number;
            
            index3=knot3.Number;
            index4=knot4.Number;
            
            dI=knot3.Value-knot1.Value;
            for i=1:numel(I2)
                I2(i).Value=I2(i).Value-dI;
            end
            %% mixing u vectors
            L1=I1(1:index1-1);
            M1=I1(index1:index2);
            R1=I1(index2+1:end);
            L2=I2(1:index3-1);
            M2=I2(index3:index4);
            R2=I2(index4+1:end);
            i=1;
            while(true)
                if L1(1).Value<L2(1).Value
                    p1{i}=L1(1);
                    if numel(L1)==1
                        status='L2';
                        break
                    end
                    L1=L1(2:end);
                    i=i+1;
                else
                    p1{i}=L2(1);
                    if numel(L2)==1
                        status='L1';
                        break
                    end
                    L2=L2(2:end);
                    i=i+1;
                end
            end
            switch status
                case 'L1'
                    for i=1:numel(L1)
                        p1{end+1}=L1(i);
                    end
                case 'L2'
                    for i=1:numel(L2)
                        p1{end+1}=L2(i);
                    end
            end
            %
            i=1;
            while(true)
                if R1(1).Value<R2(1).Value
                    p3{i}=R1(1);
                    if numel(R1)==1
                        status='R2';
                        break
                    end
                    R1=R1(2:end);
                    i=i+1;
                else
                    p3{i}=R2(1);
                    if numel(R2)==1
                        status='R1';
                        break
                    end
                    R2=R2(2:end);
                    i=i+1;
                end
            end
            switch status
                case 'R1'
                    for i=1:numel(R1)
                        p3{end+1}=R1(i);
                    end
                case 'R2'
                    for i=1:numel(R2)
                        p3{end+1}=R2(i);
                    end
            end
            %
            i=1;
            while(true)
                if M1(1).Value==M2(1).Value
                    p2{i}=M1(1);
                    if numel(M1)==1 && numel(M2)==1
                        break
                    end
                    M1=M1(2:end);
                    M2=M2(2:end);
                    i=i+1;
                elseif M1(1).Value<M2(1).Value
                    p2{i}=M1(1);
                    M1=M1(2:end);
                    i=i+1;
                else
                    p2{i}=M2(1);
                    
                    M2=M2(2:end);
                    i=i+1;
                end
            end
            %
            wholeknotU=Uknot.empty(numel(p1)+numel(p2)+numel(p3),0);
            for i=1:numel(p1)
                wholeknotU(i)=p1{i};
            end
            c=i;
            for i=c+1:c+numel(p2)
                wholeknotU(i)=p2{i-c};
            end
            c=i;
            for i=c+1:c+numel(p3)
                wholeknotU(i)=p3{i-c};
            end
            
            %% end of mixing U vectors
            %% this part we find the boundry controlpoints
            sharecp1={};
            point=point1;
            i=1;
            while(true)
                if point==point2
                    switch  C
                        case 0
                            sharecp1{i}=point;
                        case 1
                            
                            sharecp1{i,1}=point.LeftPoint;
                            sharecp1{i,2}=point;
                        case 2
                            sharecp1{i,1}=point.LeftPoint.LeftPoint;
                            sharecp1{i,2}=point.LeftPoint;
                            sharecp1{i,3}=point;
                    end
                    break
                end
                switch  C
                    case 0
                        sharecp1{i}=point;
                    case 1
                        
                        sharecp1{i,1}=point.LeftPoint;
                        sharecp1{i,2}=point;
                    case 2
                        sharecp1{i,1}=point.LeftPoint.LeftPoint;
                        sharecp1{i,2}=point.LeftPoint;
                        sharecp1{i,3}=point;
                end
                point=point.DownPoint;
                i=i+1;
            end
            %%
            sharecp2={};
            point=point3;
            i=1;
            while(true)
                if point==point4
                    switch  C
                        case 0
                            sharecp2{i}=point;
                        case 1
                            
                            sharecp2{i,2}=point.RightPoint;
                            sharecp2{i,1}=point;
                        case 2
                            sharecp2{i,3}=point.RightPoint.RightPoint;
                            sharecp2{i,2}=point.RightPoint;
                            sharecp2{i,1}=point;
                    end
                    break
                end
                switch  C
                    case 0
                        sharecp2{i}=point;
                    case 1
                        
                        sharecp2{i,2}=point.RightPoint;
                        sharecp2{i,1}=point;
                    case 2
                        sharecp2{i,3}=point.RightPoint.RightPoint;
                        sharecp2{i,2}=point.RightPoint;
                        sharecp2{i,1}=point;
                end
                point=point.DownPoint;
                i=i+1;
            end
            %% end of finding boundry controlpoints
            %% finding shared boundry uknots
            shareknot2={};
            i=1;
            point=point3;
            while(true)
                shareknot2{i}=point.CPCU;
                if point==point4
                    break
                end
                i=i+1;
                point=point.DownPoint;
            end
            
            shareknot1={};
            i=1;
            point=point1;
            while(true)
                shareknot1{i}=point.CPCU;
                if point==point2
                    break
                end
                i=i+1;
                point=point.DownPoint;
            end
            %% set new coordinates to new controlpoint
            switch continuity
                case  'C0'
                    for i=1:numel(sharecp1)
                        sharecp1{i}.X=(sharecp1{i}.X+sharecp2{i}.X)/2;
                        sharecp1{i}.Y=(sharecp1{i}.Y+sharecp2{i}.Y)/2;
                        sharecp1{i}.Z=(sharecp1{i}.Z+sharecp2{i}.Z)/2;
                        sharecp1{i}.W=(sharecp1{i}.W+sharecp2{i}.W)/2;
                    end
                case  'C1'
                    for i=1:size(sharecp1,1)
                        sharecp1{i,1}.X=(sharecp1{i,1}.X+sharecp2{i,1}.X)/2;
                        sharecp1{i,1}.Y=(sharecp1{i,1}.Y+sharecp2{i,1}.Y)/2;
                        sharecp1{i,1}.Z=(sharecp1{i,1}.Z+sharecp2{i,1}.Z)/2;
                        sharecp1{i,1}.W=(sharecp1{i,1}.W+sharecp2{i,1}.W)/2;
                        
                        sharecp1{i,2}.X=(sharecp1{i,2}.X+sharecp2{i,2}.X)/2;
                        sharecp1{i,2}.Y=(sharecp1{i,2}.Y+sharecp2{i,2}.Y)/2;
                        sharecp1{i,2}.Z=(sharecp1{i,2}.Z+sharecp2{i,2}.Z)/2;
                        sharecp1{i,2}.W=(sharecp1{i,2}.W+sharecp2{i,2}.W)/2;
                    end
                case  'C2'
                    for i=1:size(sharecp1,1)
                        sharecp1{i,1}.X=(sharecp1{i,1}.X+sharecp2{i,1}.X)/2;
                        sharecp1{i,1}.Y=(sharecp1{i,1}.Y+sharecp2{i,1}.Y)/2;
                        sharecp1{i,1}.Z=(sharecp1{i,1}.Z+sharecp2{i,1}.Z)/2;
                        sharecp1{i,1}.W=(sharecp1{i,1}.W+sharecp2{i,1}.W)/2;
                        
                        sharecp1{i,2}.X=(sharecp1{i,2}.X+sharecp2{i,2}.X)/2;
                        sharecp1{i,2}.Y=(sharecp1{i,2}.Y+sharecp2{i,2}.Y)/2;
                        sharecp1{i,2}.Z=(sharecp1{i,2}.Z+sharecp2{i,2}.Z)/2;
                        sharecp1{i,2}.W=(sharecp1{i,2}.W+sharecp2{i,2}.W)/2;
                        
                        sharecp1{i,3}.X=(sharecp1{i,3}.X+sharecp2{i,3}.X)/2;
                        sharecp1{i,3}.Y=(sharecp1{i,3}.Y+sharecp2{i,3}.Y)/2;
                        sharecp1{i,3}.Z=(sharecp1{i,3}.Z+sharecp2{i,3}.Z)/2;
                        sharecp1{i,3}.W=(sharecp1{i,3}.W+sharecp2{i,3}.W)/2;
                    end
            end
            
            %%
            switch C
                case 0
                    for i=1:numel(sharecp2)-1
                        ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i},'LastControlPoint',sharecp2{i+1}); %vertical edge of obj2
                        delete(ed);
                    end
                    
                    for i=1:numel(sharecp2)
                        ed=findobj(obj2.edge,'LastControlPoint',sharecp2{i},'Direction','horizontal');
                        delete(ed);
                    end
                    
                    ed=findobj(obj2.edge,'LastControlPoint',sharecp2{1},'Direction','vertical');
                    delete(ed);
                    
                    ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{end},'Direction','vertical');
                    delete(ed)
                    %%
                    for i=1:numel(sharecp1)
                        groupEdge=findobj(obj.edge,'FirstControlPoint',sharecp1{i},'Direction','horizontal');
                        for j=1:numel(groupEdge)
                            delete(groupEdge(j));
                        end
                    end
                    %%
                    for i=1:numel(sharecp2)
                        %
                        groupEdge=findobj(obj2.edge,'DirectionParameterCoordinate',sharecp2{i}.CPCU);
                        if ~isempty(groupEdge)
                            for j=1:numel(groupEdge)
                                groupEdge(j).DirectionParameterCoordinate=sharecp1{i}.CPCU;
                            end
                        end
                        %
                        groupEdge=findobj(obj2.edge,'FirstControlPoint',sharecp2{i},'Direction','horizontal');
                        for j=1:numel(groupEdge)
                            groupEdge(j).FirstTip=sharecp1{i}.CPCV;
                            groupEdge(j).FirstControlPoint=sharecp1{i};
                            sharecp1{i}.RightEdge=groupEdge(j);
                        end
                        %
                    end
                    %%
                    for i=1:numel(sharecp1)
                        group=findobj(obj2.controlpoint,'LeftPoint',sharecp2{i});
                        if ~isempty(group)
                            group.LeftPoint=sharecp1{i};
                            group.Vvector(2)=sharecp1{i}.CPCV;
                            group.Vvector(1)=sharecp1{i}.Vvector(2);
                            
                            sharecp1{i}.RightPoint=group;
                            sharecp1{i}.Vvector(4)=group.CPCV;
                            sharecp1{i}.Vvector(5)=group.RightPoint.CPCV;
                        end
                        group=findobj(obj2.controlpoint,'DownPoint',sharecp2{i});
                        if ~isempty(group)
                            group.DownPoint=sharecp1{i};
                        end
                        group=findobj(obj2.controlpoint,'UpPoint',sharecp2{i});
                        if ~isempty(group)
                            group.UpPoint=sharecp1{i};
                        end
                    end
                case 1
                    for i=1:size(sharecp2,1)-1
                        ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,1},'LastControlPoint',sharecp2{i+1,1}); %vertical edge of obj2
                        delete(ed);
                        ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,2},'LastControlPoint',sharecp2{i+1,2}); %vertical edge of obj2
                        delete(ed);
                    end
                    
                    for i=1:size(sharecp2,1)
                        ed=findobj(obj2.edge,'LastControlPoint',sharecp2{i,1},'Direction','horizontal');
                        delete(ed);
                        ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,1},'LastControlPoint',sharecp2{i,2});
                        delete(ed);
                    end
                    
                    ed=findobj(obj2.edge,'LastControlPoint',sharecp2{1,1},'Direction','vertical');
                    delete(ed);
                    ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{end,1},'Direction','vertical');
                    delete(ed)
                     ed=findobj(obj2.edge,'LastControlPoint',sharecp2{1,2},'Direction','vertical');
                    delete(ed);
                    ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{end,2},'Direction','vertical');
                    delete(ed)
                    %%
                    for i=1:size(sharecp1,1)
                        groupEdge=findobj(obj.edge,'FirstControlPoint',sharecp1{i,2},'Direction','horizontal');
                        delete(groupEdge);
                    end
                    
                    for i=1:size(sharecp2,1)
                        %
                        groupEdge=findobj(obj2.edge,'DirectionParameterCoordinate',sharecp2{i,1}.CPCU);
                        for j=1:numel(groupEdge)
                            groupEdge(j).DirectionParameterCoordinate=sharecp1{i,1}.CPCU;
                        end             
                    end                 
                    %%
                    for i=1:size(sharecp1,1)
                        if ~isempty(sharecp2{i,2}.RightEdge)
                            
                            sharecp1{i,2}.RightEdge=sharecp2{i,2}.RightEdge;
                            sharecp1{i,2}.RightEdge.FirstControlPoint=sharecp1{i,2};
                            sharecp1{i,2}.RightEdge.FirstTip=sharecp1{i,2}.CPCV;                        
                            
                            
                            sharecp1{i,2}.RightPoint=sharecp2{i,2}.RightPoint;
                            sharecp2{i,2}.RightPoint.LeftPoint=sharecp1{i,2};
                        end
                    end
                case 2
                    for i=1:size(sharecp2,1)
                        ed=findobj(obj2.edge,'LastControlPoint',sharecp2{i,1},'Direction','horizontal');
                        delete(ed);
                        ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,1},'LastControlPoint',sharecp2{i,2});
                        delete(ed);
                        
                        if i<size(sharecp2,1)
                            ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,1},'LastControlPoint',sharecp2{i+1,1});
                            delete(ed);
                            %
                            ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{i,2},'LastControlPoint',sharecp2{i+1,2});
                            delete(ed);
                        end
                        %
                        ed=findobj(obj.edge,'FirstControlPoint',sharecp1{i,3},'Direction','horizontal');
                        delete(ed)
                        %
                        ed=findobj(obj.edge,'FirstControlPoint',sharecp1{i,2},'LastControlPoint',sharecp1{i,3});
                        delete(ed);
                        %
                        
                        if i<size(sharecp2,1)
                            ed=findobj(obj.edge,'FirstControlPoint',sharecp1{i,3},'LastControlPoint',sharecp1{i+1,3});
                            delete(ed)
                        end
                    end
                    ed=findobj(obj2.edge,'LastControlPoint',sharecp2{1,1},'Direction','vertical');
                    delete(ed);
                    %
                    ed=findobj(obj2.edge,'LastControlPoint',sharecp2{1,2},'Direction','vertical');
                    delete(ed);
                    %
                    ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{end,1},'Direction','vertical');
                    delete(ed);
                    %
                    ed=findobj(obj2.edge,'FirstControlPoint',sharecp2{end,2},'Direction','vertical');
                    delete(ed);
                    %
                    ed=findobj(obj.edge,'LastControlPoint',sharecp1{1,3},'Direction','vertical');
                    delete(ed)
                    %
                    ed=findobj(obj.edge,'FirstControlPoint',sharecp1{end,3},'Direction','vertical');
                    delete(ed)
                    %
                    for i=1:size(sharecp2,1)
                        %
                        groupEdge=findobj(obj2.edge,'DirectionParameterCoordinate',sharecp2{i,1}.CPCU);
                        for j=1:numel(groupEdge)
                            groupEdge(j).DirectionParameterCoordinate=sharecp1{i,1}.CPCU;
                        end
                        %
                        sharecp1{i,2}.RightEdge=sharecp2{i,2}.RightEdge;
                        sharecp2{i,2}.RightEdge.FirstControlPoint=sharecp1{i,2};
                        sharecp2{i,2}.RightEdge.FirstTip=sharecp1{i,2}.CPCV;
                        sharecp1{i,2}.RightPoint=sharecp2{i,3};
                        sharecp2{i,3}.LeftPoint=sharecp1{i,2};
                        %
                    end
                    %%
            end
            %% this part is for merging whole side of surfaces
            for i=1:numel(obj2.controlpoint)
                for k=1:5
                    if obj2.controlpoint(i).Uvector(k)==obj2.U(1)
                        obj2.controlpoint(i).Uvector(k)=obj.U(1);
                    elseif  obj2.controlpoint(i).Uvector(k)==obj2.U(2)
                        obj2.controlpoint(i).Uvector(k)=obj.U(2);
                    elseif obj2.controlpoint(i).Uvector(k)==obj2.U(end-1)
                        obj2.controlpoint(i).Uvector(k)=obj.U(end-1);
                    elseif obj2.controlpoint(i).Uvector(k)==obj2.U(end)
                        obj2.controlpoint(i).Uvector(k)=obj.U(end);
                    end
                end
            end
            i=1;
            while(true)
                if ~isvalid(obj2.edge(i))
                    obj2.edge(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(obj2.edge)
                    break
                end
            end
            
            for i=1:numel(obj2.edge)
                %
                if obj2.edge(i).FirstTip== obj2.U(1)
                    obj2.edge(i).FirstTip=sharecp1{1}.Uvector(1);
                elseif obj2.edge(i).LastTip== obj2.U(1)
                    obj2.edge(i).LastTip=sharecp1{1}.Uvector(1);
                elseif obj2.edge(i).DirectionParameterCoordinate== obj2.U(1)
                    obj2.edge(i).DirectionParameterCoordinate=sharecp1{1}.Uvector(1);
                end
                %
                if obj2.edge(i).FirstTip== obj2.U(2)
                    obj2.edge(i).FirstTip=sharecp1{1}.Uvector(2);
                elseif obj2.edge(i).LastTip== obj2.U(2)
                    obj2.edge(i).LastTip=sharecp1{1}.Uvector(2);
                elseif obj2.edge(i).DirectionParameterCoordinate== obj2.U(2)
                    obj2.edge(i).DirectionParameterCoordinate=sharecp1{1}.uvector(2);
                end
                %
                if obj2.edge(i).FirstTip== obj2.U(end-1)
                    obj2.edge(i).FirstTip=sharecp1{end}.Uvector(4);
                elseif obj2.edge(i).LastTip== obj2.U(end-1)
                    obj2.edge(i).LastTip=sharecp1{end}.Uvector(4);
                elseif obj2.edge(i).DirectionParameterCoordinate== obj2.U(end-1)
                    obj2.edge(i).DirectionParameterCoordinate=sharecp1{end}.Uvector(4);
                end
                %
                if obj2.edge(i).FirstTip== obj2.U(end)
                    obj2.edge(i).FirstTip=sharecp1{end}.Uvector(end);
                elseif obj2.edge(i).LastTip== obj2.U(end)
                    obj2.edge(i).LastTip=sharecp1{end}.Uvector(end);
                elseif obj2.edge(i).DirectionParameterCoordinate== obj2.U(end)
                    obj2.edge(i).DirectionParameterCoordinate=sharecp1{end}.Uvector(end);
                end
            end
            %%
            i=1;
            while(true)
                if ~isvalid(obj.controlpoint(i))
                    obj.controlpoint(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(obj.controlpoint)
                    break
                end
            end
            %%
            i=1;
            while(true)
                if ~isvalid(obj.edge(i))
                    obj.edge(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(obj.edge)
                    break
                end
            end
            %%
            i=1;
            while(true)
                if ~isvalid(obj2.controlpoint(i))
                    obj2.controlpoint(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(obj2.controlpoint)
                    break
                end
            end
            %%
            i=1;
            while(true)
                if ~isvalid(obj2.edge(i))
                    obj2.edge(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(obj2.edge)
                    break
                end
            end
            %%
            for i=1:numel(shareknot1)
                for j=1:numel(obj2.controlpoint)
                    for k=1:5
                        if obj2.controlpoint(j).Uvector(k)==shareknot2{i}
                            obj2.controlpoint(j).Uvector(k)=shareknot1{i};
                        end
                    end
                end
                %%
                for j=1:numel(obj2.edge)
                    if isa(obj2.edge(j).FirstTip,'Uknot')
                        if obj2.edge(j).FirstTip==shareknot2{i}
                            obj2.edge(j).FirstTip=shareknot1{i};
                        elseif obj2.edge(j).LastTip==shareknot2{i}
                            obj2.edge(j).LastTip=shareknot1{i};
                        end
                    end
                end
                
                delete(shareknot2{i})
            end
            %%
            i=1;
            while(true)
                if ~isvalid(wholeknotU(i))
                    wholeknotU(i)=[];
                    i=i-1;
                end
                i=i+1;
                if i>numel(wholeknotU)
                    break
                end
            end
            
            %%
            wholepoint=ControlPoint.empty(numel(obj.controlpoint)+numel(obj2.controlpoint),0);
            for i=1:numel(obj.controlpoint)+numel(obj2.controlpoint)
                if i<=numel(obj.controlpoint)
                    wholepoint(i)=obj.controlpoint(i);
                else
                    wholepoint(i)=obj2.controlpoint(i-numel(obj.controlpoint));
                end
            end
            %%
            switch C
                case 0
                    for i=1:numel(sharecp2)
                        delete(sharecp2{i});
                    end
                case 1
                    for i=1:numel(sharecp2)
                        delete(sharecp2{i});
                    end
                case 2
                    for i=1:size(sharecp2,1)
                        delete(sharecp2{i,1})
                        delete(sharecp2{i,2})
                        delete(sharecp1{i,3})
                    end
            end
            %%
            wholeedge=Edge.empty(numel(obj.edge)+numel(obj2.edge),0);
            wholeedge(1:numel(obj.edge))=obj.edge;
            wholeedge(1+numel(obj.edge):numel(obj.edge)+numel(obj2.edge))=obj2.edge;
            %%
            newobj=T_Surface;
            newobj.controlpoint=wholepoint;
            newobj.U=wholeknotU;
            newobj.V=wholeknotV;
            newobj.edge=wholeedge;
            
            %%editing newobj.U , newobj.V
            intervalU=[];
            for i=1:numel(newobj.U)-1
                intervalU(i)=newobj.U(i+1).Value-newobj.U(i).Value;
            end
            
            intervalU=intervalU/sum(intervalU,2);
            newobj.U(1).Value=0;
            for i=2:numel(newobj.U)
                newobj.U(i).Value=newobj.U(i-1).Value+intervalU(i-1);
            end
            
            intervalV=[];
            for i=1:numel(newobj.V)-1
                intervalV(i)=newobj.V(i+1).Value-newobj.V(i).Value;
            end
            intervalV=intervalV/sum(intervalV,2);
            newobj.V(1).Value=0;
            for i=2:numel(newobj.V)
                newobj.V(i).Value=newobj.V(i-1).Value+intervalV(i-1);
            end
            %%
            %to delete exra Vknots
            switch C
                case 0
                    delete(obj.V(end-1:end));
                    delete(obj2.V(1:3));
                case 1
                    delete(obj.V(end-1:end));
                    delete(obj2.V(1:4));
                case 2
                    delete(obj.V(end-2:end))
                    delete(obj2.V(1:4))
            end
            %%
            obj.controlpoint=newobj.controlpoint;
            obj.edge=newobj.edge;
            obj.U=newobj.U;
            obj.V=newobj.V;
            %%
            newobj=obj;
            i=1;
            while(true)
                if ~isvalid(newobj.controlpoint(i)) || isempty(newobj.controlpoint(i))
                    newobj.controlpoint(i)=[];
                    i=i-1;
                end
                i=i+1;
                if numel(newobj.controlpoint)<i
                    break
                end
            end
            %%
            i=1;
            while(true)
                if ~isvalid(newobj.edge(i)) || isempty(newobj.edge(i))
                    newobj.edge(i)=[];
                    i=i-1;
                end
                if ~isvalid(newobj.edge(i).DirectionParameterCoordinate) || isempty(newobj.edge(i).DirectionParameterCoordinate)
                    newobj.edge(i)=[];
                    i=i-1;
                end
                i=i+1;
                if numel(newobj.edge)<i
                    break
                end
            end
            for i=1:numel(newobj.controlpoint)
                newobj.controlpoint(i).Number=i;
            end
            for i=1:numel(newobj.edge)
                newobj.edge(i).Number=i;
            end
            %
            delete(obj2.U(1:2));
            delete(obj2.U(end-1:end));
            %%
            i=1;
            while(true)
                if ~isvalid(newobj.U(i)) || isempty(newobj.U(i))
                    newobj.U(i)=[];
                    i=i-1;
                end
                i=i+1;
                if numel(newobj.U)<i
                    break
                end
            end
            for i=1:numel(newobj.U)
                newobj.U(i).Number=i;
            end
            %%
            i=1;
            while(true)
                if ~isvalid(newobj.V(i)) || isempty(newobj.V(i))
                    newobj.V(i)=[];
                    i=i-1;
                end
                i=i+1;
                if numel(newobj.V)<i
                    break
                end
            end
            for i=1:numel(newobj.V)
                newobj.V(i).Number=i;
            end
            %%
            %%backknot and nextknot
            for i=1:numel(obj.V)
                if i<numel(obj.V)
                    obj.V(i).NextKnot=obj.V(i+1);
                end
                if i>1
                    obj.V(i).BackKnot=obj.V(i-1);
                end
            end
            obj.V(1).BackKnot=[];
            obj.V(end).NextKnot=[];
            %%backknot and nextknot
            for i=1:numel(obj.U)
                if i<numel(obj.U)
                    obj.U(i).NextKnot=obj.U(i+1);
                end
                if i>1
                    obj.U(i).BackKnot=obj.U(i-1);
                end
            end
            obj.U(1).BackKnot=[];
            obj.U(end).NextKnot=[];
            %%%
            %for insurence
            for i=1:numel(newobj.controlpoint)
                [u,v]=newobj.CheckDomain(newobj.controlpoint(i).Number);
                newobj.controlpoint(i).Uvector=u;
                newobj.controlpoint(i).Vvector=v;
            end
        end
        %%
        function obj=ReverseKnotVector(obj,direction)% this function is to reverse the selected knot vector of
            % T-spline surface
            % the inputs are T-spline object and the direction of knot
            % vector('horizontal' means v vector and 'vertical' means u
            % vector)
            switch direction
                case 'U'
                    for i=1:numel(obj.U)
                        obj.U(i).Value=1-obj.U(i).Value;
                    end
                    %%sorting U
                    u=obj.U;
                    for i=1:numel(obj.U)
                        obj.U(i)=u(end-i+1);
                    end
                    %%numbering
                    for i=1:numel(obj.U)
                        obj.U(i).Number=i;
                    end
                    %%backknot and nextknot
                    for i=1:numel(obj.U)
                        if i<numel(obj.U)
                            obj.U(i).NextKnot=obj.U(i+1);
                        end
                        if i>1
                            obj.U(i).BackKnot=obj.U(i-1);
                        end
                    end
                    %%
                    for i=1:numel(obj.controlpoint)
                        u1=obj.controlpoint(i).Uvector(1);
                        u2=obj.controlpoint(i).Uvector(2);
                        u4=obj.controlpoint(i).Uvector(4);
                        u5=obj.controlpoint(i).Uvector(5);
                        
                        obj.controlpoint(i).Uvector(1)=u5;
                        obj.controlpoint(i).Uvector(2)=u4;
                        obj.controlpoint(i).Uvector(4)=u2;
                        obj.controlpoint(i).Uvector(5)=u1;
                    end
                    for i=1:numel(obj.controlpoint)
                        c=obj.controlpoint(i).UpPoint;
                        cc=obj.controlpoint(i).DownPoint;
                        
                        obj.controlpoint(i).UpPoint=cc;
                        obj.controlpoint(i).DownPoint=c;
                        %
                        c=obj.controlpoint(i).UpEdge;
                        cc=obj.controlpoint(i).DownEdge;
                        
                        obj.controlpoint(i).UpEdge=cc;
                        obj.controlpoint(i).DownEdge=c;
                        
                    end
                    %%
                    groupEdge=findobj(obj.edge,'Direction','vertical');
                    for i=1:numel(groupEdge)
                        first=groupEdge(i).FirstTip;
                        last=groupEdge(i).LastTip;
                        
                        groupEdge(i).FirstTip=last;
                        groupEdge(i).LastTip=first;
                        %
                        first=groupEdge(i).FirstControlPoint;
                        last=groupEdge(i).LastControlPoint;
                        
                        groupEdge(i).FirstControlPoint=last;
                        groupEdge(i).LastControlPoint=first;
                    end
                    
                case 'V'
                    for i=1:numel(obj.V)
                        obj.V(i).Value=1-obj.V(i).Value;
                    end
                    %%sorting V
                    v=obj.V;
                    for i=1:numel(obj.V)
                        obj.V(i)=v(end-i+1);
                    end
                    %%numbering
                    for i=1:numel(obj.V)
                        obj.V(i).Number=i;
                    end
                    %%backknot and nextknot
                    for i=1:numel(obj.V)
                        if i<numel(obj.V)
                            obj.V(i).NextKnot=obj.V(i+1);
                        end
                        if i>1
                            obj.V(i).BackKnot=obj.V(i-1);
                        end
                    end
                    
                    
                    %%
                    for i=1:numel(obj.controlpoint)
                        c=obj.controlpoint(i).LeftPoint;
                        cc=obj.controlpoint(i).RightPoint;
                        
                        obj.controlpoint(i).LeftPoint=cc;
                        obj.controlpoint(i).RightPoint=c;
                        %
                        c=obj.controlpoint(i).LeftEdge;
                        cc=obj.controlpoint(i).RightEdge;
                        
                        obj.controlpoint(i).LeftEdge=cc;
                        obj.controlpoint(i).RightEdge=c;
                        
                        v1=obj.controlpoint(i).Vvector(1);
                        v2=obj.controlpoint(i).Vvector(2);
                        v4=obj.controlpoint(i).Vvector(4);
                        v5=obj.controlpoint(i).Vvector(5);
                        
                        obj.controlpoint(i).Vvector(1)=v5;
                        obj.controlpoint(i).Vvector(2)=v4;
                        obj.controlpoint(i).Vvector(4)=v2;
                        obj.controlpoint(i).Vvector(5)=v1;
                    end
                    %%
                    groupEdge=findobj(obj.edge,'Direction','horizontal');
                    for i=1:numel(groupEdge)
                        first=groupEdge(i).FirstTip;
                        last=groupEdge(i).LastTip;
                        
                        groupEdge(i).FirstTip=last;
                        groupEdge(i).LastTip=first;
                        %
                        first=groupEdge(i).FirstControlPoint;
                        last=groupEdge(i).LastControlPoint;
                        
                        groupEdge(i).FirstControlPoint=last;
                        groupEdge(i).LastControlPoint=first;
                    end
            end
        end
        %%
        function obj=SwitchVectors(obj) %this function is to switch or exchange the knot vector togather. this function is mostly used in merging function
            newU=Uknot.empty(numel(obj.V),0);
            for i=1:numel(obj.V)
                newU(i)=Uknot;
                newU(i).Value=obj.V(i).Value;
                newU(i).Number=obj.V(i).Number;
            end
            %
            for i=1:numel(newU)
                if i<numel(newU)
                    newU(i).NextKnot=newU(i+1);
                end
                if i>1
                    newU(i).BackKnot=newU(i-1);
                end
            end
            
            newV=Vknot.empty(numel(obj.U),0);
            for i=1:numel(obj.U)
                newV(i)=Vknot;
                newV(i).Value=obj.U(i).Value;
                newV(i).Number=obj.U(i).Number;
            end
            %
            for i=1:numel(newV)
                if i<numel(newV)
                    newV(i).NextKnot=newV(i+1);
                end
                if i>1
                    newV(i).BackKnot=newV(i-1);
                end
            end
            %
            reserveknot=cell(numel(obj.controlpoint),2);
            
            for i=1:size(reserveknot,1)
                reserveknot{i,1}=cell(1,5);
                reserveknot{i,2}=cell(1,5);
            end
            %
            for i=1:numel(newU)
                for j=1:numel(obj.controlpoint)
                    for k=1:5
                        if obj.controlpoint(j).Vvector(k).Number==i
                            reserveknot{j,1}{1,k}=newU(i);
                        end
                    end
                end
            end
            
            for i=1:numel(newV)
                for j=1:numel(obj.controlpoint)
                    for k=1:5
                        if obj.controlpoint(j).Uvector(k).Number==i
                            reserveknot{j,2}{1,k}=newV(i);
                        end
                    end
                end
            end
            
            for i=1:numel(obj.controlpoint)
                u=Uknot.empty(5,0);
                v=Vknot.empty(5,0);
                for j=1:5
                    u(j)=reserveknot{i,1}{1,j};
                    v(j)=reserveknot{i,2}{1,j};
                end
                obj.controlpoint(i).Uvector=u;
                obj.controlpoint(i).Vvector=v;
            end
            %
            for i=1:numel(obj.edge)
                if strcmp(obj.edge(i).Direction,'horizontal')
                    u1=findobj(newU,'Number',obj.edge(i).FirstTip.Number);
                    obj.edge(i).FirstTip=u1;
                    
                    u2=findobj(newU,'Number',obj.edge(i).LastTip.Number);
                    obj.edge(i).LastTip=u2;
                    
                    v=findobj(newV,'Number',obj.edge(i).DirectionParameterCoordinate.Number);
                    obj.edge(i).DirectionParameterCoordinate=v;
                    
                    obj.edge(i).Direction='vertical';
                else
                    v1=findobj(newV,'Number',obj.edge(i).FirstTip.Number);
                    obj.edge(i).FirstTip=v1;
                    
                    v2=findobj(newV,'Number',obj.edge(i).LastTip.Number);
                    obj.edge(i).LastTip=v2;
                    
                    u=findobj(newU,'Number',obj.edge(i).DirectionParameterCoordinate.Number);
                    obj.edge(i).DirectionParameterCoordinate=u;
                    
                    obj.edge(i).Direction='horizontal';
                end
            end
            %
            delete(obj.V(:));
            delete(obj.U(:));
            
            obj.V=newV;
            obj.U=newU;
            
            %
            for i=1:numel(obj.controlpoint)
                c=obj.controlpoint(i).RightPoint;
                cc=obj.controlpoint(i).DownPoint;
                %
                obj.controlpoint(i).RightPoint=cc;
                obj.controlpoint(i).DownPoint=c;
                %%%
                k=obj.controlpoint(i).UpPoint;
                kk=obj.controlpoint(i).LeftPoint;
                %
                obj.controlpoint(i).UpPoint=kk;
                obj.controlpoint(i).LeftPoint=k;
                %
                %%%%%%%%%
                c=obj.controlpoint(i).RightEdge;
                cc=obj.controlpoint(i).DownEdge;
                %
                obj.controlpoint(i).RightEdge=cc;
                obj.controlpoint(i).DownEdge=c;
                %
                k=obj.controlpoint(i).UpEdge;
                kk=obj.controlpoint(i).LeftEdge;
                %
                obj.controlpoint(i).UpEdge=kk;
                obj.controlpoint(i).LeftEdge=k;
            end
        end
        %%
        function plotmeshTp(obj) %this function plot the topolgy of T-spline surface      
            %the axis of this plot is index of knots in to vu and v
            %parameter
            figure(2)
            format short
            digits(2)
            hold on
            
            %%%
            for i=1:numel(obj.edge)
                try
                    %                 if ~isempty(obj.edge(i).FirstControlPoint) && ~isempty(obj.edge(i).LastControlPoint)
                    if strcmp(obj.edge(i).Direction,'horizontal')
                        line([obj.edge(i).FirstTip.Number,obj.edge(i).LastTip.Number],...
                            [obj.edge(i).DirectionParameterCoordinate.Number,obj.edge(i).DirectionParameterCoordinate.Number],'color','k');
                        text((obj.edge(i).FirstTip.Number+obj.edge(i).LastTip.Number)/2....
                            , obj.edge(i).DirectionParameterCoordinate.Number+0.1....
                            , [ num2str(obj.edge(i).Number),'-',char(vpa(obj.edge(i).LastTip.Value-obj.edge(i).FirstTip.Value))],'color','k','HorizontalAlignment','center','fontsize',10)
                    else
                        line([obj.edge(i).DirectionParameterCoordinate.Number,obj.edge(i).DirectionParameterCoordinate.Number],.....
                            [obj.edge(i).FirstTip.Number,obj.edge(i).LastTip.Number],'color','k');
                        text(obj.edge(i).DirectionParameterCoordinate.Number+0.05,(obj.edge(i).FirstTip.Number+obj.edge(i).LastTip.Number)/2,....
                            [ num2str(obj.edge(i).Number),'-',char(vpa(obj.edge(i).LastTip.Value-obj.edge(i).FirstTip.Value))],'color','k','HorizontalAlignment','center','rotation',-90,'fontsize',10);
                    end
                    %                 end
                catch
                    warning(['edge i:',num2str(i)])
                    continue
                end
            end
            xlabel('knot number of V parameter')
            ylabel('knot number of U parameter')
            
            for i=1:numel(obj.controlpoint)
                try
                    plot(obj.controlpoint(i).CPCV.Number,obj.controlpoint(i).CPCU.Number,'o','markerfacecolor','y','markerfacecolor','w','markeredgecolor','w')
                    text(obj.controlpoint(i).CPCV.Number,obj.controlpoint(i).CPCU.Number,num2str(obj.controlpoint(i).Number),'color',[0.5 0  1],'fontsize',12,'fontweight','bold','HorizontalAlignment','center')
                catch
                    warning(['controlpoint i:',num2str(i)])
                    continue
                end
            end
        end
        %%
        function [coordinateX,coordinateY,coordinateZ]=Calculate(obj,Uparameter,Vparameter) %the purpose of this function is to calculate of a point on the T-spline surface. 
            %the inputs are T-spline object, u and v parameter of the
            %point.
            p=3;% degree
            u=Uparameter;
            v=Vparameter;
            
            for i=1:numel(obj.controlpoint)
                if u>=obj.controlpoint(i).Uvector(1).Value && u<obj.controlpoint(i).Uvector(end).Value
                    if v>=obj.controlpoint(i).Vvector(1).Value && v<obj.controlpoint(i).Vvector(end).Value 
                        [N,Q]=T_Surface.basisfunction(obj.controlpoint(i),u,v,p);
                        X=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).X*obj.controlpoint(i).W;
                        Y=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).Y*obj.controlpoint(i).W;
                        Z=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).Z*obj.controlpoint(i).W;
                        Weight=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).W;
                        if ~exist('x','var')
                            x=X;
                            y=Y;
                            z=Z;
                            w=Weight;
                        else
                            x=x+X;
                            y=y+Y;
                            z=z+Z;
                            w=w+Weight;
                        end
                    end
                end
            end
            
            %%for insurance
            
            if exist('x','var')
                %%
                coordinateX=x/w;
                coordinateY=y/w;
                coordinateZ=z/w;
                
                %%
                clear x y z
            else
                coordinateX=[];
                coordinateY=[];
                coordinateZ=[];
            end
        end
        %%
        function [dux,duy,duz,dvx,dvy,dvz]=Derivate1(obj,Uparameter,Vparameter)%the derivation of surface is base on bspline theory but not
            %nurbs. it means it must be developed by nurbs   
            
            p=3;% degree
            u=Uparameter;
            v=Vparameter;
            
            for i=1:numel(obj.controlpoint)
                if u>=obj.controlpoint(i).Uvector(1).Value && u<obj.controlpoint(i).Uvector(end).Value
                    if v>=obj.controlpoint(i).Vvector(1).Value && v<obj.controlpoint(i).Vvector(end).Value
                        [N,Q]=T_Surface.basisfunction(obj.controlpoint(i),u,v,p);
                        %%
                        DuL=p*N(1,p)/(obj.controlpoint(i).Uvector(4).Value-obj.controlpoint(i).Uvector(1).Value);
                        if isnan(DuL)==1 || isinf(DuL)==1
                            DuL=0;
                        end
                        DuR=p*N(2,p)/(obj.controlpoint(i).Uvector(5).Value-obj.controlpoint(i).Uvector(2).Value);
                        if isnan(DuR)==1 || isinf(DuR)==1
                            DuR=0;
                        end
                        Dux=(DuL-DuR)*Q(1,p+1)*obj.controlpoint(i).X*obj.controlpoint(i).W;
                        Duy=(DuL-DuR)*Q(1,p+1)*obj.controlpoint(i).Y*obj.controlpoint(i).W;
                        Duz=(DuL-DuR)*Q(1,p+1)*obj.controlpoint(i).Z*obj.controlpoint(i).W;
                        %%
                        DvL=p*Q(1,p)/(obj.controlpoint(i).Vvector(4).Value-obj.controlpoint(i).Vvector(1).Value);
                        if isnan(DvL)==1 || isinf(DvL)==1
                            DvL=0;
                        end
                        DvR=p*Q(2,p)/(obj.controlpoint(i).Vvector(5).Value-obj.controlpoint(i).Vvector(2).Value);
                        if isnan(DvR)==1 || isinf(DvR)==1
                            DvR=0;
                        end
                        Dvx=(DvL-DvR)*N(1,p+1)*obj.controlpoint(i).X*obj.controlpoint(i).W;
                        Dvy=(DvL-DvR)*N(1,p+1)*obj.controlpoint(i).Y*obj.controlpoint(i).W;
                        Dvz=(DvL-DvR)*N(1,p+1)*obj.controlpoint(i).Z*obj.controlpoint(i).W;
                        
                        Weight=N(1,p+1)*Q(1,p+1)*obj.controlpoint(i).W;
                        if ~exist('dux','var')
                            dux=Dux;
                            duy=Duy;
                            duz=Duz;
                            
                            dvx=Dvx;
                            dvy=Dvy;
                            dvz=Dvz;
                            w=Weight;
                        else
                            dux=dux+Dux;
                            duy=duy+Duy;
                            duz=duz+Duz;
                            
                            dvx=dvx+Dvx;
                            dvy=dvy+Dvy;
                            dvz=dvz+Dvz;
                            
                            w=w+Weight;
                        end
                    end
                end
            end
            
            if ~exist('dux','var')
                dux=[];
                duy=[];
                duz=[];
                
                dvx=[];
                dvy=[];
                dvz=[];
            end
        end
        %%
        function [UvectorDetail,VvectorDetail]=UVDetail(obj) % the purpose of the this function is to show the knot vectors and their details in table in command window
            for i=1:numel(obj.U)
                if ~isvalid(obj.U(i))
                    UNum{i}='---';
                    UValue{i}='---';
                    BackU{i}='---';
                    NextU{i}='---';
                else
                    UNum{i}=num2str(obj.U(i).Number);
                    UValue{i}=num2str(obj.U(i).Value);
                    if isempty(obj.U(i).BackKnot)
                        BackU{i}='EMPTY';
                    else
                        BackU{i}=num2str(obj.U(i).BackKnot.Number);
                    end
                    
                    if isempty(obj.U(i).NextKnot)
                        NextU{i}='EMPTY';
                    else
                        NextU{i}=num2str(obj.U(i).NextKnot.Number);
                    end
                end
            end
            UNum=UNum';UValue=UValue';BackU=BackU';NextU=NextU';
            UvectorDetail=table(UNum,UValue,BackU,NextU);
            %
            for i=1:numel(obj.V)
                if ~isvalid(obj.V(i))
                    VNum{i}='---';
                    VValue{i}='---';
                    BackV{i}='---';
                    NextV{i}='---';
                else
                    VNum{i}=num2str(obj.V(i).Number);
                    VValue{i}=num2str(obj.V(i).Value);
                    if isempty(obj.V(i).BackKnot)
                        BackV{i}='EMPTY';
                    else
                        BackV{i}=num2str(obj.V(i).BackKnot.Number);
                    end
                    
                    if isempty(obj.V(i).NextKnot)
                        NextV{i}='EMPTY';
                    else
                        NextV{i}=num2str(obj.V(i).NextKnot.Number);
                    end
                end
            end
            VNum=VNum';VValue=VValue';BackV=BackV';NextV=NextV';
            VvectorDetail=table(VNum,VValue,BackV,NextV);
            %
        end
        %%
        function obj=ReArrange(obj) %make controlpoints and vectors in order by number
            for i=1:numel(obj.controlpoint)
                obj.controlpoint(i).Number=i;
            end
            for i=1:numel(obj.edge)
                obj.edge(i).Number=i;
            end
            for i=1:numel(obj.U)
                obj.U(i).Number=i;
            end
            for i=1:numel(obj.V)
                obj.V(i).Number=i;
            end
        end
        %%
        function [udomain,vdomain]=CheckDomain(obj,numberCP)
            vdomain=Vknot.empty(5,0);
            udomain=Uknot.empty(5,0);
            %finding controlpoint
            cp=findobj(obj.controlpoint,'Number',numberCP);
            %finding Vdomain
            %%%the area of point
                vdomain(3)=cp.CPCV;
                c=0;
                for i=cp.CPCV.Number+1:obj.V(end).Number
                    
                    if c==2
                        break
                    end
                    if i>obj.V(end).Number-3
                        c=c+1;
                        vdomain(3+c)=obj.V(i);
                        continue
                    end
                    p=findobj(obj.controlpoint,'CPCU',cp.CPCU,'CPCV',obj.V(i));
                    if ~isempty(p)
                        c=c+1;
                        vdomain(3+c)=obj.V(i);
                        continue
                    end
                    ed=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(i));
                    if ~isempty(ed)
                        for j=1:numel(ed)
                            if ed(j).FirstTip.Number<cp.CPCU.Number && ed(j).LastTip.Number>cp.CPCU.Number
                                c=c+1;
                                vdomain(3+c)=obj.V(i);
                                break
                            end
                        end
                    end
                end
                %
                c=0;
                for i=cp.CPCV.Number-1:-1:1
                    if c==2
                        break
                    end
                     if i<3
                        c=c+1;
                        vdomain(3-c)=obj.V(i);
                        continue
                    end
                    p=findobj(obj.controlpoint,'CPCU',cp.CPCU,'CPCV',obj.V(i));
                    if ~isempty(p)
                        c=c+1;
                        vdomain(3-c)=obj.V(i);
                        continue
                    end
                    
                    ed=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(i));
                    if ~isempty(ed)
                        for j=1:numel(ed)
                            if ed(j).FirstTip.Number<=cp.CPCU.Number && ed(j).LastTip.Number>=cp.CPCU.Number
                                c=c+1;
                                vdomain(3-c)=obj.V(i);
                                break
                            end
                        end
                    end
                end
            %%%Udomain checking
            udomain(3)=cp.CPCU;
                c=0;
                for i=cp.CPCU.Number+1:obj.U(end).Number
                    
                    if c==2
                        break
                    end
                     if i>obj.U(end).Number-3
                        c=c+1;
                        udomain(3+c)=obj.U(i);
                        continue
                    end
                    p=findobj(obj.controlpoint,'CPCV',cp.CPCV,'CPCU',obj.U(i));
                    if ~isempty(p)
                        c=c+1;
                        udomain(3+c)=obj.U(i);
                        continue
                    end
                    
                    ed=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(i));
                    if ~isempty(ed)
                        for j=1:numel(ed)
                            if ed(j).FirstTip.Number<=cp.CPCV.Number && ed(j).LastTip.Number>=cp.CPCV.Number
                                c=c+1;
                                udomain(3+c)=obj.U(i);
                                break
                            end
                        end
                    end 
                end
                %
                c=0;
                for i=cp.CPCU.Number-1:-1:1
                    if c==2
                        break
                    end
                     if i<3
                        c=c+1;
                        udomain(3-c)=obj.U(i);
                        continue
                    end
                    p=findobj(obj.controlpoint,'CPCV',cp.CPCV,'CPCU',obj.U(i));
                    if ~isempty(p)
                        c=c+1;
                        udomain(3-c)=obj.U(i);
                        continue
                    end
                    
                    ed=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(i));
                    if ~isempty(ed)
                        for j=1:numel(ed)
                            if ed(j).FirstTip.Number<=cp.CPCV.Number && ed(j).LastTip.Number>=cp.CPCV.Number
                                c=c+1;
                                udomain(3-c)=obj.U(i);
                                break
                            end
                        end
                    end
                end
            %%%
            vdomain(3)=cp.CPCV;
            udomain(3)=cp.CPCU;
            
        end%%%checking controlpoint Domain by mesh
        %%
        function obj=InsertRow(obj,knot) %to make a row or collumn of control points on the control grid
            % the inputs are T-spline object and a knot
            % it's important that the knot must be existed. if you want to
            % add a knot you must use inserting point function to build or
            % add a knot.
            if isa(knot,'Vknot')
                for i=3:numel(obj.U)-2
                    knotU=obj.U(i);
                    cp=findobj(obj.controlpoint,'CPCU',knotU,'CPCV',knot);
                    if isempty(cp)
                        ed=findobj(obj.edge,'DirectionParameterCoordinate',knotU);
                        for j=1:numel(ed)
                            if knot.Number<ed(j).LastTip.Number  &&  knot.Number>ed(j).FirstTip.Number
                                obj=InsertPoint(obj,ed(j).Number,knot);
                                break
                            end
                        end
                    end
                end
                
            else
                for i=3:numel(obj.V)-2
                    knotV=obj.V(i);
                    cp=findobj(obj.controlpoint,'CPCU',knot,'CPCV',knotV);
                    if isempty(cp)
                        ed=findobj(obj.edge,'DirectionParameterCoordinate',knotV);
                        for j=1:numel(ed)
                            if knot.Number<ed(j).LastTip.Number  &&  knot.Number>ed(j).FirstTip.Number
                                obj=InsertPoint(obj,ed(j).Number,knot);
                                break
                            end
                        end
                    end
                end
            end
        end
        %%
        function [varargout]=plotDomainCP(obj,cp_num)
            d=0.1;
            
            Right=line([obj.controlpoint(cp_num).Vvector(5).Number+d,obj.controlpoint(cp_num).Vvector(5).Number+d]....
                ,[obj.controlpoint(cp_num).Uvector(1).Number-d,obj.controlpoint(cp_num).Uvector(5).Number+d],'linewidth',3,'color','r');
            
            Left=line([obj.controlpoint(cp_num).Vvector(1).Number-d,obj.controlpoint(cp_num).Vvector(1).Number-d]....
                ,[obj.controlpoint(cp_num).Uvector(1).Number-d,obj.controlpoint(cp_num).Uvector(5).Number+d],'linewidth',3,'color','r');
            
            Up=line([obj.controlpoint(cp_num).Vvector(1).Number-d,obj.controlpoint(cp_num).Vvector(5).Number+d]....
                ,[obj.controlpoint(cp_num).Uvector(1).Number-d,obj.controlpoint(cp_num).Uvector(1).Number-d],'linewidth',3,'color','r');
            
            Down=line([obj.controlpoint(cp_num).Vvector(1).Number-d,obj.controlpoint(cp_num).Vvector(5).Number+d]....
                ,[obj.controlpoint(cp_num).Uvector(5).Number+d,obj.controlpoint(cp_num).Uvector(5).Number+d],'linewidth',3,'color','r');
             P=plot(obj.controlpoint(cp_num).CPCV.Number,obj.controlpoint(cp_num).CPCU.Number,'o','markerfacecolor','b','markersize',10);
            if nargout==1
                varargout{1}={};
                varargout{1}{1}=Right;
                varargout{1}{2}=Down;
                varargout{1}{3}=Up;
                varargout{1}{4}=Left;
                varargout{1}{5}=P;
            end
            
        end %to show domain of control point in Topology(you must plot topology then use this function)
        %%
        function plotmeshp(obj) % the purposse of this function is to plot control grid in parametric coordinate
            figure(3)
            format short
            digits(2)
            hold on
            for i=1:numel(obj.controlpoint)
                try
                    plot(obj.controlpoint(i).CPCV.Value,obj.controlpoint(i).CPCU.Value,'o','markerfacecolor','y','markerfacecolor','k','markeredgecolor','k')
%                     text(obj.controlpoint(i).CPCV.Value,obj.controlpoint(i).CPCU.Value,num2str(obj.controlpoint(i).Number),'color',[0.5 0  1],'fontsize',15,'fontweight','bold','HorizontalAlignment','center')
                catch
                    warning(['controlpoint i:',num2str(i)])
                    continue
                end
            end
            %%%
            for i=1:numel(obj.edge)
                try
                    if ~isempty(obj.edge(i).FirstControlPoint) && ~isempty(obj.edge(i).LastControlPoint)
                        if strcmp(obj.edge(i).Direction,'horizontal')
                            line([obj.edge(i).FirstTip.Value,obj.edge(i).LastTip.Value],...
                                [obj.edge(i).DirectionParameterCoordinate.Value,obj.edge(i).DirectionParameterCoordinate.Value],'color','k');
                            
                        else
                            line([obj.edge(i).DirectionParameterCoordinate.Value,obj.edge(i).DirectionParameterCoordinate.Value],.....
                                [obj.edge(i).FirstTip.Value,obj.edge(i).LastTip.Value],'color','k');
                            
                        end
                    end
                catch
                    warning(['edge i:',num2str(i)])
                    continue
                end
            end
            xlabel('v value')
            ylabel('u value')
            axis equal
            %             axis equal
        end
        %%
        function [A,B]=Coef(obj,Point,u,v,Knot,Direction)%the purpose of this function is to calculate the linear combination of Basis Function of a selected control point
            %% inputs %%
            % obj: T-spline surface
            % Point: seleceted control point
            % u: primary u vector
            % v: primary v vector
            % Knot: additional knot in u or v vector
            % Direction: 'horizontal' or 'vertical'
            if strcmp(Direction,'horizontal')
                Vector=v;
            else
                Vector=u;
            end
            % calculating the coeficients   c1 , c2
            if   Vector(1).Number < Knot.Number &&  Vector(2).Number > Knot.Number
                c1=(Knot.Value-Vector(1).Value)/(Vector(4).Value-Vector(1).Value);
                c2=1;
                count=1;
            elseif Vector(2).Number < Knot.Number &&  Vector(3).Number > Knot.Number
                c1=(Knot.Value-Vector(1).Value)/(Vector(4).Value-Vector(1).Value);
                c2=(Vector(5).Value-Knot.Value)/(Vector(5).Value-Vector(2).Value);
                count=2;
            elseif Vector(3).Number < Knot.Number &&  Vector(4).Number > Knot.Number
                c1=(Knot.Value-Vector(1).Value)/(Vector(4).Value-Vector(1).Value);
                c2=(Vector(5).Value-Knot.Value)/(Vector(5).Value-Vector(2).Value);
                count=3;
            elseif Vector(4).Number < Knot.Number &&  Vector(5 ).Number > Knot.Number
                c1=1;
                c2=(Vector(5).Value-Knot.Value)/(Vector(5).Value-Vector(2).Value);
                count=4;
            else
                error('not in vector')
            end
%%
            wholeV={};
            for i=1:6
                if i<=count
                    wholeV{end+1}=Vector(i);
                elseif i==count+1
                    wholeV{end+1}=Knot;
                else
                    wholeV{end+1}=Vector(i-1);
                end
            end
            %%
            vector1=[wholeV{1},wholeV{2},wholeV{3},wholeV{4},wholeV{5}];
            vector2=[wholeV{2},wholeV{3},wholeV{4},wholeV{5},wholeV{6}];
            %%
            if strcmp(Direction,'horizontal')             
                cp1=findobj(obj.controlpoint,'CPCU',Point.CPCU,'CPCV',vector1(3));
                cp2=findobj(obj.controlpoint,'CPCU',Point.CPCU,'CPCV',vector2(3));
                Vvector1=vector1; Vvector2=vector2;Uvector1=u;Uvector2=u;
            else
                cp1=findobj(obj.controlpoint,'CPCU',vector1(3),'CPCV',Point.CPCV);
                cp2=findobj(obj.controlpoint,'CPCU',vector2(3),'CPCV',Point.CPCV);
                Uvector1=vector1; Uvector2=vector2;Vvector1=v;Vvector2=v;
            end
            %%
            A={}; B={}; 
            A{1}=Uvector1;A{2}=Vvector1;A{3}=c1;A{4}=cp1;
            B{1}=Uvector2;B{2}=Vvector2;B{3}=c2;B{4}=cp2;
            
        end
        %%
        function obj=Remesh(obj,Point) %this function is to remesh the control grid by adding a new point it's mainly used to show former algorithm to add control points
            %%
            for i=Point.CPCV.Number-1:-1:1
                EDs=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(i));
                status='unfinished';
                for j=1:numel(EDs)
                    if EDs(j).FirstTip.Number<=Point.CPCU.Number && EDs(j).LastTip.Number>=Point.CPCU.Number
                        CP=findobj(obj.controlpoint,'CPCU',Point.CPCU,'CPCV',obj.V(i));
                        ed=findobj(obj.edge,'FirstControlPoint',CP,'Direction','horizontal');
                        if ~isempty(CP) && isempty(ed)
                            ed=Edge('horizontal');
                            ed.FirstTip=CP.CPCV;    ed.LastTip=Point.CPCV;     ed.FirstControlPoint=CP;     ed.LastControlPoint=Point;
                            ed.DirectionParameterCoordinate=Point.CPCU; ed.Direction='horizontal';
                            obj.edge=[obj.edge,ed];
                            ed.Number=numel(obj.edge);
                            
                            Point.LeftEdge=ed; Point.LeftPoint=CP; 
                            CP.RightEdge=ed; CP.RightPoint=Point;
                            status='finished';
                            break
                        end
                    end
                end
                if strcmp(status,'finished')
                    break
                end
                if ~isempty(EDs) && strcmp(status,'unfinished')
                    break
                end
            end
            %%
            for i=Point.CPCV.Number+1:obj.V(end).Number
                EDs=findobj(obj.edge,'DirectionParameterCoordinate',obj.V(i));
                status='unfinished';
                for j=1:numel(EDs)
                    if EDs(j).FirstTip.Number<=Point.CPCU.Number && EDs(j).LastTip.Number>=Point.CPCU.Number
                        CP=findobj(obj.controlpoint,'CPCU',Point.CPCU,'CPCV',obj.V(i));
                        ed=findobj(obj.edge,'LastControlPoint',CP,'Direction','horizontal');
                        if ~isempty(CP) && isempty(ed)
                            ed=Edge('horizontal');
                            ed.FirstTip=Point.CPCV;    ed.LastTip=CP.CPCV;     ed.FirstControlPoint=Point;     ed.LastControlPoint=CP;
                            ed.DirectionParameterCoordinate=Point.CPCU; ed.Direction='horizontal';
                            
                            obj.edge=[obj.edge,ed];
                            ed.Number=numel(obj.edge);
                            
                            Point.RightEdge=ed; Point.RightPoint=CP; 
                            CP.LeftEdge=ed; CP.LeftPoint=Point;
                            status='finished';
                            break
                        end
                    end
                end
                if strcmp(status,'finished')
                    break
                end
                if ~isempty(EDs) && strcmp(status,'unfinished')
                    break
                end
            end
            %%
            for i=Point.CPCU.Number-1:-1:1
                EDs=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(i));
                status='unfinished';
                for j=1:numel(EDs)
                    if EDs(j).FirstTip.Number<=Point.CPCV.Number && EDs(j).LastTip.Number>=Point.CPCV.Number
                        CP=findobj(obj.controlpoint,'CPCU',obj.U(i),'CPCV',Point.CPCV);
                        ed=findobj(obj.edge,'FirstControlPoint',CP,'Direction','vertical');
                        if ~isempty(CP) && isempty(ed)
                            ed=Edge('vertical');
                            ed.FirstTip=CP.CPCU;    ed.LastTip=Point.CPCU;     ed.FirstControlPoint=CP;     ed.LastControlPoint=Point;
                            ed.DirectionParameterCoordinate=Point.CPCV; ed.Direction='vertical';
                            
                            obj.edge=[obj.edge,ed];
                            ed.Number=numel(obj.edge);
                            
                            Point.UpEdge=ed; Point.UpPoint=CP; 
                            CP.DownEdge=ed; CP.DownPoint=Point;
                            status='finished';
                            break
                        end
                    end
                end
                if strcmp(status,'finished')
                    break
                end
                if ~isempty(EDs) && strcmp(status,'unfinished')
                    break
                end
            end
            %%
            for i=Point.CPCU.Number+1:obj.U(end).Number
                EDs=findobj(obj.edge,'DirectionParameterCoordinate',obj.U(i));
                status='unfinished';
                for j=1:numel(EDs)
                    if EDs(j).FirstTip.Number<=Point.CPCV.Number && EDs(j).LastTip.Number>=Point.CPCV.Number
                        CP=findobj(obj.controlpoint,'CPCU',obj.U(i),'CPCV',Point.CPCV);
                        ed=findobj(obj.edge,'LastControlPoint',CP,'Direction','vertical');
                        if ~isempty(CP) && isempty(ed)
                            ed=Edge('vertical');
                            ed.FirstTip=Point.CPCU;    ed.LastTip=CP.CPCU;     ed.FirstControlPoint=Point;     ed.LastControlPoint=CP;
                            ed.DirectionParameterCoordinate=Point.CPCV; ed.Direction='vertical';
                            
                            obj.edge=[obj.edge,ed];
                            ed.Number=numel(obj.edge);
                            
                            Point.DownEdge=ed; Point.DownPoint=CP; 
                            CP.UpEdge=ed; CP.UpPoint=Point;
                            status='finished';
                            break
                        end
                    end
                end
                if strcmp(status,'finished')
                    break
                end
                if ~isempty(EDs) && strcmp(status,'unfinished')
                    break
                end
            end
            %%
            
            eds = findobj(obj.edge,'DirectionParameterCoordinate',Point.CPCV);
            if ~isempty(eds)
                for i=1:numel(eds)
                    if Point.CPCU.Number < eds(i).LastTip.Number && Point.CPCU.Number > eds(i).FirstTip.Number
                        newedge = Edge('vertical');
                        newedge.FirstTip = Point.CPCU; newedge.LastTip = eds(i).LastTip;
                        newedge.FirstControlPoint = Point; newedge.LastControlPoint = eds(i).LastControlPoint;
                        newedge.DirectionParameterCoordinate = Point.CPCV;
                        %
                        eds(i).LastTip = Point.CPCU; eds(i).LastControlPoint = Point;
                        %
                        eds(i).FirstControlPoint.DownPoint = eds(i).FirstControlPoint;
                        newedge.LastControlPoint.UpPoint = newedge.FirstControlPoint;
                        %
                        Point.UpPoint = eds(i).FirstControlPoint;
                        eds(i).FirstControlPoint.DownPoint = Point;
                        
                        Point.DownPoint = newedge.LastControlPoint;
                        newedge.LastControlPoint.UpPoint = Point;
                        %
                        obj.edge = [obj.edge,newedge];
                        newedge.Number = numel(obj.edge);
                    end
                end
            end
            %
            eds = findobj(obj.edge,'DirectionParameterCoordinate',Point.CPCU);
            if ~isempty(eds)
                for i=1:numel(eds)
                    if Point.CPCV.Number < eds(i).LastTip.Number && Point.CPCV.Number > eds(i).FirstTip.Number
                        newedge = Edge('horizontal');
                        newedge.FirstTip = Point.CPCV; newedge.LastTip = eds(i).LastTip;
                        newedge.FirstControlPoint = Point; newedge.LastControlPoint = eds(i).LastControlPoint;
                        newedge.DirectionParameterCoordinate = Point.CPCU;
                        %
                        eds(i).LastTip = Point.CPCV; eds(i).LastControlPoint = Point;
                        %
                        eds(i).FirstControlPoint.RightPoint = eds(i).FirstControlPoint;
                        newedge.LastControlPoint.LeftPoint = newedge.FirstControlPoint;
                        %
                        Point.LeftPoint = eds(i).FirstControlPoint;
                        eds(i).FirstControlPoint.RightPoint = Point;
                        
                        Point.RightPoint = newedge.LastControlPoint;
                        newedge.LastControlPoint.LeftPoint = Point;
                        %
                        obj.edge = [obj.edge,newedge];
                        newedge.Number = numel(obj.edge);
                    end
                end
            end
            
            
        end
        %%
        function obj=C0Curves(obj)     %this function calculate the c0 continuty curves on surface  

            
            %%C0 curves in V direction
             C0vertical=[];  % to show parametric range of c0 curves on surface 1th collum is v parameter an 2th and 3ht collums are Umin and Umax 
            for i=3:numel(obj.U)-3
                
                gp=findobj(obj.edge,'Direction','vertical');
                groupEdge=Edge.empty;
                for j=1:numel(gp)
                    if gp(j).FirstTip.Number <= obj.U(i).Number && gp(j).LastTip.Number > obj.U(i).Number
                        groupEdge(end+1)=gp(j);
                    end
                end
                
                if numel(groupEdge)>2
                    j=1;
                    while(true)
                        if groupEdge(j).DirectionParameterCoordinate.Number > groupEdge(j+1).DirectionParameterCoordinate.Number
                            s=groupEdge(j);
                            groupEdge(j)=groupEdge(j+1);
                            groupEdge(j+1)=s;
                            j=j-2;
                        end
                        j=j+1;
                        if j==numel(groupEdge)
                            break
                        end
                    end
                    
                    for k=1:numel(groupEdge)
                        if k+2<=numel(groupEdge) && groupEdge(k).DirectionParameterCoordinate.Value-groupEdge(k+1).DirectionParameterCoordinate.Value==0  .....
                                && groupEdge(k+1).DirectionParameterCoordinate.Value-groupEdge(k+2).DirectionParameterCoordinate.Value==0
                            C0vertical(end+1,1)=groupEdge(k+1).DirectionParameterCoordinate.Value;
                            C0vertical(end,2)=obj.U(i).Value;
                            C0vertical(end,3)=obj.U(i+1).Value;
                            
                        end
                    end
                end
            end
            %%%
            %%C0 curves in U direction
             C0horizontal=[];  % to show parametric range of c0 curves on surface 1th collum is v parameter an 2th and 3ht collums are Umin and Umax 
            for i=3:numel(obj.V)-3
                                              
                gp=findobj(obj.edge,'Direction','horizontal');
                groupEdge=Edge.empty;
                for j=1:numel(gp)
                    if gp(j).FirstTip.Number <= obj.V(i).Number && gp(j).LastTip.Number > obj.V(i).Number
                        groupEdge(end+1)=gp(j);
                    end
                end
                
                if numel(groupEdge)>2
                    j=1;
                    while(true)
                        if groupEdge(j).DirectionParameterCoordinate.Number > groupEdge(j+1).DirectionParameterCoordinate.Number
                            s=groupEdge(j);
                            groupEdge(j)=groupEdge(j+1);
                            groupEdge(j+1)=s;
                            j=j-2;
                        end
                        j=j+1;
                        if j==numel(groupEdge)
                            break
                        end
                    end
                    
                    for k=1:numel(groupEdge)
                        if k+2<=numel(groupEdge) && groupEdge(k).DirectionParameterCoordinate.Value-groupEdge(k+1).DirectionParameterCoordinate.Value==0  .....
                                && groupEdge(k+1).DirectionParameterCoordinate.Value-groupEdge(k+2).DirectionParameterCoordinate.Value==0
                            C0horizontal(end+1,1)=groupEdge(k+1).DirectionParameterCoordinate.Value;
                            C0horizontal(end,2)=obj.V(i).Value;
                            C0horizontal(end,3)=obj.V(i+1).Value;
                            
                        end
                    end
                end  
            end
            %%sorting of  C0horizontal and C0vertical based on u and v
            
            if size(C0vertical,1)>1 
                i=1;
                while(true)
                    if C0vertical(i,1) > C0vertical(i+1,1)
                        s=C0vertical(i+1,:);
                        C0vertical(i+1,:)=C0vertical(i,:);
                        C0vertical(i,:)=s;
                        i=i-1;
                    end
                    i=i+1;
                    if i==size(C0vertical,1)
                        break
                    end
                end
            end
            %%
             if size(C0horizontal,1)>1 
                i=1;
                while(true)
                    if C0horizontal(i,1) > C0horizontal(i+1,1)
                        s=C0horizontal(i+1,:);
                        C0horizontal(i+1,:)=C0horizontal(i,:);
                        C0horizontal(i,:)=s;
                        i=i-1;
                    end
                    i=i+1;
                    if i==size(C0horizontal,1)
                        break
                    end
                end
             end
            %%end of sorting
            curveU=[]; curveV=[];
            %%
            if ~isempty(C0vertical)
            curveV(end+1,1)=C0vertical(1,1); curveV(end,2)=C0vertical(1,2); curveV(end,3)=C0vertical(1,3);
            for i=2:size(C0vertical,1)
                if curveV(end,1)==C0vertical(i,1)
                    if curveV(end,3)==C0vertical(i,2)
                        curveV(end,3)=C0vertical(i,3);
                    else
                        curveV(end+1,1)=C0vertical(i,1);
                        curveV(end,2)=C0vertical(i,2);
                        curveV(end,3)=C0vertical(i,3);
                    end
                else
                    curveV(end+1,1)=C0vertical(i,1);
                    curveV(end,2)=C0vertical(i,2);
                    curveV(end,3)=C0vertical(i,3);
                end
            end
            end
            %
            if ~isempty(C0horizontal)
            curveU(end+1,1)=C0horizontal(1,1); curveU(end,2)=C0horizontal(1,2); curveU(end,3)=C0horizontal(1,3);
            for i=2:size(C0horizontal,1)
                if curveU(end,1)==C0horizontal(i,1)
                    if curveU(end,3)==C0horizontal(i,2)
                        curveU(end,3)=C0horizontal(i,3);
                    else
                        curveU(end+1,1)=C0horizontal(i,1);
                        curveU(end,2)=C0horizontal(i,2);
                        curveU(end,3)=C0horizontal(i,3);
                    end
                else
                    curveU(end+1,1)=C0horizontal(i,1);
                    curveU(end,2)=C0horizontal(i,2);
                    curveU(end,3)=C0horizontal(i,3);
                end
            end
            end
            %
            curveData={};
            
            if ~isempty(curveV)
                for i=1:size(curveV,1)
                    u=linspace(curveV(i,2),curveV(i,3),100);
                    points=zeros(100,3);
                    for j=1:100
                        [x,y,z]=Calculate(obj,u(j),curveV(i,1));
                        if ~isempty(x)
                            points(j,1)=x;
                            points(j,2)=y;
                            points(j,3)=z;
                        end
                    end
                    curveData{end+1}=points;
                end
            end
            if ~isempty(curveU)
                for i=1:size(curveU,1)
                    v=linspace(curveU(i,2),curveU(i,3),100);
                    points=zeros(100,3);
                    for j=1:100
                        [x,y,z]=Calculate(obj,curveU(i,1),v(j));
                        if ~isempty(x)
                            points(j,1)=x;
                            points(j,2)=y;
                            points(j,3)=z;
                        end
                    end
                    curveData{end+1}=points;
                end
            end
            obj.C0_points=curveData;
            obj.C0_curveU=curveU;
            obj.C0_curveV=curveV;
        end 
        %%
        function a=Partition_Of_Unity(obj,u,v) %this function is to check partion of unity in specific parameters.
            a=0;
            for i=1:numel(obj.controlpoint)
                if u>=obj.controlpoint(i).Uvector(1).Value && u<obj.controlpoint(i).Uvector(end).Value
                    if  v>=obj.controlpoint(i).Vvector(1).Value && v<obj.controlpoint(i).Vvector(end).Value
                        [N,Q]=T_Surface.basisfunction(obj.controlpoint(i),u,v,3);
                        a=a+N(1,end)*Q(1,end);
                    end
                end
            end
        end
        %%
        function [varargout]=Neighbor(obj,cp_num)
            d=0.1;
            point=obj.controlpoint(cp_num);
            %%
            if ~isempty(point.RightEdge)
                Right_Edge=line([point.RightEdge.FirstTip.Number,point.RightEdge.LastTip.Number]....
                    ,[point.RightEdge.DirectionParameterCoordinate.Number,point.RightEdge.DirectionParameterCoordinate.Number],'linewidth',3,'color','r');
            end
            if ~isempty(point.LeftEdge)
                Left_Edge=line([point.LeftEdge.FirstTip.Number,point.LeftEdge.LastTip.Number]....
                    ,[point.LeftEdge.DirectionParameterCoordinate.Number,point.LeftEdge.DirectionParameterCoordinate.Number],'linewidth',3,'color','r');
            end
            if ~isempty(point.UpEdge)
                Up_Edge=line([point.UpEdge.DirectionParameterCoordinate.Number,point.UpEdge.DirectionParameterCoordinate.Number]....
                    ,[point.UpEdge.FirstTip.Number,point.UpEdge.LastTip.Number],'linewidth',3,'color','r');
            end
            if ~isempty(point.DownEdge)
                Down_Edge=line([point.DownEdge.DirectionParameterCoordinate.Number,point.DownEdge.DirectionParameterCoordinate.Number]....
                    ,[point.DownEdge.FirstTip.Number,point.DownEdge.LastTip.Number],'linewidth',3,'color','r');
            end
            P=plot(point.CPCV.Number,point.CPCU.Number,'o','markerfacecolor','b','markersize',10);
            %%
            if ~isempty(point.RightPoint)
                R=plot(point.RightPoint.CPCV.Number,point.RightPoint.CPCU.Number,'ro','markerfacecolor','r','markersize',10);
            end
            if ~isempty(point.LeftPoint)
                L=plot(point.LeftPoint.CPCV.Number,point.LeftPoint.CPCU.Number,'ro','markerfacecolor','r','markersize',10);
            end
            if ~isempty(point.UpPoint)
                U=plot(point.UpPoint.CPCV.Number,point.UpPoint.CPCU.Number,'ro','markerfacecolor','r','markersize',10);
            end
            if ~isempty(point.DownPoint)
                D=plot(point.DownPoint.CPCV.Number,point.DownPoint.CPCU.Number,'ro','markerfacecolor','r','markersize',10);
            end
            %%
            varargout{1}={};
            if nargout==1
                
                if exist('Right_Edge')
                    varargout{1}{end+1}=Right_Edge;
                end
                if exist('Left_Edge')
                    varargout{1}{end+1}=Left_Edge;
                end
                if exist('Up_Edge')
                    varargout{1}{end+1}=Up_Edge;
                end
                if exist('Down_Edge')
                    varargout{1}{end+1}=Down_Edge;
                end
                %%
                if exist('R')
                    varargout{1}{end+1}=R;
                end
                if exist('L')
                    varargout{1}{end+1}=L;
                end
                if exist('U')
                    varargout{1}{end+1}=U;
                end
                if exist('D')
                    varargout{1}{end+1}=D;
                end
                 varargout{1}{end+1}=P;
               end
        end       %this function to show the neighboring edges and control points of a specific control point     
        %%       
    end
    methods(Static=true)
        function [N,Q]=basisfunction(controlpoint,u,v,p)
            NN=zeros(p+1);
            QQ=zeros(p+1);
            for L=1:p+1
                if u>=controlpoint.Uvector(L).Value && u<controlpoint.Uvector(L+1).Value
                    i=L;
                end
                if v>=controlpoint.Vvector(L).Value && v<controlpoint.Vvector(L+1).Value
                    j=L;
                end
            end
            NN(i,1)=1;
            QQ(j,1)=1;
            for q=1:p
                for j=1:p+1-q
                    leftN=(u-controlpoint.Uvector(j).Value)/(controlpoint.Uvector(j+q).Value-controlpoint.Uvector(j).Value);
                    rightN=(controlpoint.Uvector(j+q+1).Value-u)/(controlpoint.Uvector(j+q+1).Value-controlpoint.Uvector(j+1).Value);
                    if isinf(leftN) || isnan(leftN)
                        leftN=0;
                    end
                    if isinf(rightN) || isnan(rightN)
                        rightN=0;
                    end
                    leftQ=(v-controlpoint.Vvector(j).Value)/(controlpoint.Vvector(j+q).Value-controlpoint.Vvector(j).Value);
                    rightQ=(controlpoint.Vvector(j+q+1).Value-v)/(controlpoint.Vvector(j+q+1).Value-controlpoint.Vvector(j+1).Value);
                    if isinf(leftQ) || isnan(leftQ)
                        leftQ=0;
                    end
                    if isinf(rightQ) || isnan(rightQ)
                        rightQ=0;
                    end
                    NN(j,q+1)=leftN.*NN(j,q)+rightN.*NN(j+1,q);
                    QQ(j,q+1)=leftQ.*QQ(j,q)+rightQ.*QQ(j+1,q);
                end
            end
          if i ==1 && u == 0
              NN(1,end) = 1;
          elseif i == 1 && u == 1 
              NN(1,end) = 1;
          end
           if j == 1 && v == 0
              QQ(1,end) = 1;
          elseif j == 1 && v == 1 
              QQ(1,end) = 1;
          end
            N=NN;
            Q=QQ;
        end%his fnction is to calculate basis function of specific control point. it's mostly used for point calculation, point derivation ad other functions 
    end
end