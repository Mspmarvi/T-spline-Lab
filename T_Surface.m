classdef T_Surface<handle
    
    properties%(all of them are array of objects)
        controlpoint %vector of ControlPoint objects
        edge %vector of Edge objects
        U   %vector of Uknot objects
        V   %vector of Vknot objects
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
                        text((obj.edge(i).FirstControlPoint.X+obj.edge(i).LastControlPoint.X)/2,....
                            (obj.edge(i).FirstControlPoint.Y+obj.edge(i).LastControlPoint.Y)/2,...
                            (obj.edge(i).FirstControlPoint.Z+obj.edge(i).LastControlPoint.Z)/2,....
                            num2str(obj.edge(i).Number),'Color',[1 0 0],'Fontsize',11);
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
                    plot3(obj.controlpoint(i).X,obj.controlpoint(i).Y,obj.controlpoint(i).Z,'o','markersize',6,'markerfacecolor','b','markeredgecolor','b');
                end
            catch
                error(['i= ',num2str(i)]);
            end
            %%
            try
                %this part is to draw control point in plot
                for i=1:numel(obj.controlpoint)
                    text(obj.controlpoint(i).X,obj.controlpoint(i).Y,obj.controlpoint(i).Z,num2str(obj.controlpoint(i).Number),'Fontsize',12,'color','k')
                end
            catch
                error(['i= ',num2str(i)]);
            end
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