function [] = PlotI3DterationHistory(nor,NodeCoord,IC,nelx,nely,nelz,...
    scale,PhiVals,order)
%% Description
% This function is used for plot. users can change views 
% or axis order for getting desired view
%% Inputs
% nor = Number of rectangle
% NodeCoord = Node coordinates
% IC variable is used for prohibitting duplicate indexing, because 
    % nodes with the same coordinates in different rectangle cubes must 
    % have unique global indexes
% nelx,nely and nelz are nummber of elements for each rectangle cube at 
    % each direction
% PhiVals saves Phi values at nodes for each iteration
% scale is a number which multiples number of nodes for plotting result
% order represents order of axis for plot, for example order = [3,1,2]
%   replace z and x axis, also order must be equal to [1,2,3] or [2,3,1] or
%   [3,1,2]
%%
arguments
    nor,NodeCoord,IC,nelx,nely,nelz,scale,PhiVals
    order = [2,3,1];
end
if ismember(order,[1,2,3;2,3,1;3,1,2],"rows")
%==========================================================================
    % InputCoord is used for saving inputted coordiantes
    % QueryCoord is (nor X 3) cell which saves query coordinates of scaled
    % coordinates
    % p1 and p2 are patches used for plot
    QueryCoord = cell(nor,3);
    InputCoord = cell(nor,3);
    p1 = cell(nor,1);
    p2 = cell(nor,1);
    for i=1:nor
    
        y = NodeCoord(:,2);
        x = NodeCoord(:,1);
        z = NodeCoord(:,3);
        
        InputCoord{i,1} = reshape(y(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
        InputCoord{i,2} = reshape(x(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
        InputCoord{i,3} = reshape(z(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
    
        y2 = linspace(min(InputCoord{i,1}(:)), max(InputCoord{i,1}(:)),...
            scale*(nely(i)+1));
        x2 = linspace(min(InputCoord{i,2}(:)), max(InputCoord{i,2}(:)),...
            scale*(nelx(i)+1));
        z2 = linspace(min(InputCoord{i,3}(:)), max(InputCoord{i,3}(:)),...
            scale*(nelz(i)+1));
    
        [QueryCoord{i,1}, QueryCoord{i,2}, QueryCoord{i,3}] = ...
            meshgrid(y2, x2, z2);
    
        p1{i} = patch('FaceColor','r','EdgeColor','none');
        p2{i} = patch('FaceColor','r','EdgeColor','none');
    
    end
    % ax = axes();
    camlight 'right' 
    view([36 7]) 
    axis 'equal'
    axis off
    for j=1:numel(PhiVals)
        for i=1:nor
            phi = reshape(PhiVals{j}(IC{i}), (nelx(i)+1), (nely(i)+1), (nelz(i)+1));
            phiq = interpn(InputCoord{i,2},...
                           InputCoord{i,1},...
                           InputCoord{i,3},...
                           phi,...
                           QueryCoord{i,2},...
                           QueryCoord{i,1},...
                           QueryCoord{i,3});
            [ff,vf] = isosurface(QueryCoord{i,order(1)},...
                                 QueryCoord{i,order(2)},...
                                 QueryCoord{i,order(3)},...
                                 phiq,0);
            [fc,vc,~] = isocaps(QueryCoord{i,order(1)},...
                QueryCoord{i,order(2)},...
                QueryCoord{i,order(3)},...
                phiq,0);
            set(p1{i}, 'Faces',ff,'Vertices',vf)
            set(p2{i}, 'Faces',fc,'Vertices',vc)
        end
        drawnow
    end
else
    ME = MException('MyComponent:noSuchVariable', ...
        'order must be equal to [1,2,3] or [2,3,1] or [3,1,2]');
    throw(ME)
end
end

