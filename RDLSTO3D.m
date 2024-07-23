%**************************************************************************
% MeshGeneration
[ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,...
    nely,nelz] = MeshGeneration(GP.nor,GP.lx,GP.ly,GP.lz,GP.els,...
    GP.X0,GP.Y0,GP.Z0,GP.BF);
%**************************************************************************
% Boundary conditions implementation
[F,Dis,FixedNodes,ukdis,kdis] = BoundaryConditionsImplementation(ndof,...
    GP.EBC,GP.NBC,NodeCoord,GP.els);
%**************************************************************************
% Indexes which is used in stiffness assembly
iK = reshape(kron(EDofMat,ones(24,1))',24*24*noe,1);
jK = reshape(kron(EDofMat,ones(1,24))',24*24*noe,1);
%**************************************************************************
% Element stiffness calculation
[Ke,Vole,Bmat,Emat] = BrickElementStiffnessCalculation(GP.ngp,GP.els,MP.v);
%**************************************************************************
% Calculation of T1 and T2 for each element
[T1e,T2e] = RD_T1T2Generator(GP.ngp,GP.els);
% T1 and T2 assembly
iT = reshape(kron(ElNodes,ones(8,1))',8*8*noe,1);
jT = reshape(kron(ElNodes,ones(1,8))',8*8*noe,1);
st1 = reshape(T1e(:)*ones(1,noe),8*8*noe,1);
st2 = reshape(T2e(:)*ones(1,noe),8*8*noe,1);
T1 = sparse(iT,jT,st1);
T2 = sparse(iT,jT,st2);
%**************************************************************************
% Vonmises stress calculation at vertices
% [Vonmises] = Vonmises_stress(noe,Bmat_node,Dis,EDofMat,Emat,...
%     Emax,NodeRep,ElNodes);
%**************************************************************************
%% Optimization loop
%**************************************************************************
% Parameters
% opc = Optimization counter
opc = 1; 
% [r,s,t] = Grid on each element for Phi interpolation
[r,s,t] = meshgrid(linspace(-1,1,20)); 
% Compliance is objective function
Compliance = zeros(OP.NOI,1);
% VC is volume constraint
VC = zeros(OP.NOI,1);
% if OP.PlotResults
%     figure
% end
PhiVals = cell(OP.NOI,1);
%**************************************************************************
while opc<=OP.NOI
    % Volume fraction calculation
    % tmpPhi = Interpolated Phi values on element for volume fraction
    % calculation
    tmpPhi = 0.125*((1-r(:)).*(1-s(:)).*(1-t(:))*Phi(ElNodes(:,1))'+...
                    (1+r(:)).*(1-s(:)).*(1-t(:))*Phi(ElNodes(:,2))'+...
                    (1+r(:)).*(1+s(:)).*(1-t(:))*Phi(ElNodes(:,3))'+...
                    (1-r(:)).*(1+s(:)).*(1-t(:))*Phi(ElNodes(:,4))'+...
                    (1-r(:)).*(1-s(:)).*(1+t(:))*Phi(ElNodes(:,5))'+...
                    (1+r(:)).*(1-s(:)).*(1+t(:))*Phi(ElNodes(:,6))'+...
                    (1+r(:)).*(1+s(:)).*(1+t(:))*Phi(ElNodes(:,7))'+...
                    (1-r(:)).*(1+s(:)).*(1+t(:))*Phi(ElNodes(:,7))');
    % vlfe = volume fraction of each element which can be obtained by
    % sum(tmpPhi>=0)/numel(r or s or t)
    vlfe = sum(tmpPhi>=0)/numel(r);
    VC(opc) = sum(vlfe)/noe;
%**************************************************************************
    % ErModel=Ersatz material model and stiffness calculation
    ErModel = vlfe*MP.Emax+(1-vlfe)*MP.Emin;
%**************************************************************************
    % Assembling of stiffness matrix
    sK = MP.Emax*reshape(Ke(:)*ErModel,24*24*noe,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2; % for assuring symmetry
%**************************************************************************
    % Solving equlibrium equations
    Dis(ukdis) = K(ukdis,ukdis)\(F(ukdis)-K(ukdis,kdis)*Dis(kdis));
%**************************************************************************
    % ElemComp = Elemet compliance
    ElemComp=sum(0.5*(Ke*Dis(EDofMat)').*(Dis(EDofMat)'.*ErModel));
    Compliance(opc) = sum(ElemComp);
%**************************************************************************
    % Check convegence
    if opc>OP.nRelax && abs(VC(opc)-OP.VC0)/OP.VC0<1e-5 && ...
            all(abs(Compliance(opc)-Compliance(opc-9:opc-1))...
            /Compliance(opc)<1e-5)
        break;
    end
%**************************************************************************
    % Comp_phi = objective function sensitivity wrt Phi
    Comp_phi = sparse(ElNodes,ones(noe,8),0.125*ElemComp'.*ones(noe,8));
%**************************************************************************
    % Lambda calcultion in augmented Lagrngian
    if opc<OP.nRelax
        Lambda=OP.Mu*(VC(opc)-VC(1)+(VC(1)-OP.VC0)*opc/OP.nRelax);
    else
        Lambda =Lambda+OP.Gamma*(VC(opc)-OP.VC0);
        OP.Gamma = min(OP.Gamma+OP.dGamma,OP.maxGamma);
    end
%**************************************************************************
    % V = Boundary velocity in LS method
    V = -Comp_phi/mean(abs(Comp_phi))+Lambda;
%**************************************************************************
    % Updating scheme in RD level set method
    T = (T1/OP.delta_T+OP.Tho*T2);
    Y = (T1*(Phi/OP.delta_T-V));
    Phi=T\Y;
    Phi = min(max(Phi,-1),1);
    Phi(FixedNodes) = 1;
%**************************************************************************
    % PhiVals saves Phi values at nodes for each iteration
    PhiVals{opc} = Phi;
%**************************************************************************
    [VC(opc),opc]
    opc = opc+1;
end
%% Mesh generation function
function [ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,...
    nely,nelz] = MeshGeneration(nor,lx,ly,lz,els,X0,Y0,Z0,BF)
%==========================================================================
    % Inputs
    % nor = Number of rectangles
    % lx = array which has x_length of rectangles
    % ly = array which has y_length of rectangles
    % lz = array which has z_length of rectangles
    % els = Element size at each direction
%==========================================================================
    % Preallocation arrays for the 
    % nummber of elements for each rectangle at each direction
    nelx = zeros(nor,1); % nelx(i) = int64(lx(i)/els); 
    nely = zeros(nor,1); % nely(i) = int64(ly(i)/els); 
    nelz = zeros(nor,1); % nelz(i) = int64(lz(i)/els);
%==========================================================================
    % coordinates{i} = [x,y,z] where x,y and z are node coordinates for
    % i'th rectangle
    coordinates = cell(nor,1);
%==========================================================================
    % NodeCoord = [coordinates{1};...,coordinates{nor}] by removing
    % duplicate nodes
    NodeCoord = [];
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nelx(i) = int64(lx(i)/els); 
        nely(i) = int64(ly(i)/els); 
        nelz(i) = int64(lz(i)/els);
        nonx = nelx(i)+1; 
        nony = nely(i)+1; 
        nonz = nelz(i)+1;
%==========================================================================
        % X0 = X-coordinate of left-bottom-back for each rectangle
        % Y0 = Y-coordinate of left-bottom-back for each rectangle
        % Z0 = Z-coordinate of left-bottom-back for each rectangle
        % nonx, nony and nonz are number of nodes for i'th rectangle for
        % x-y-z direction, respectively which can be calculated by
        x_node = linspace(X0(i),X0(i)+lx(i),nonx);
        y_node = linspace(Y0(i),Y0(i)+ly(i),nony);
        z_node = linspace(Z0(i),Z0(i)+lz(i),nonz);
%==========================================================================
        x = repmat(x_node',nony*nonz,1);
        y = repmat(kron(y_node',ones(nonx,1)),nonz,1);
        z = kron(z_node',ones(nonx*nony,1));
        coordinates{i} = [x,y,z];
        NodeCoord = [NodeCoord;coordinates{i}];
    end
%==========================================================================
    % removing duplicate nodes
    NodeCoord = unique(NodeCoord,'stable','rows');
%==========================================================================
    % Phi specifies if material exist at nodes. size(Phi)=number of rows in
    % NodeCoord
    % Phi>=0 => material exist
    % Phi = 0=>boundary representation
    % Phi<0 => no material
    Phi= 0.5*ones(size(NodeCoord,1),1);
%==========================================================================
    % Loop over BF for boundary specification(where Phi=0)
        for j=1:size(BF,1)
            ind = NodeCoord(:,1)>=BF(j,1) & NodeCoord(:,1)<=BF(j,2)...
                & NodeCoord(:,2)>=BF(j,3) & NodeCoord(:,2)<=BF(j,4)...
                & NodeCoord(:,3)>=BF(j,5) & NodeCoord(:,3)<=BF(j,6);
            Phi(ind)=0;
        end 
%==========================================================================
    % ElNodes = global indexes for nodes of each element 
    ElNodes = [];
%==========================================================================
    % IC variable for prohibitting duplicate indexing is used, because 
    % nodes with the same coordinates in different rectangles must have 
    % unique global indexes     
    IC = cell(nor,1);
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nodes = (1:size(coordinates{i},1))';
        [~,ic,iC] = intersect(coordinates{i},NodeCoord,'stable','rows');
        nodes(ic) = iC;IC{i}=iC;
        nonx = nelx(i)+1; nony = nely(i)+1;
        a = repmat([0,1,nelx(i)+2,nelx(i)+1],nelx(i),1);
        b = repmat(a,nely(i),1)+kron((0:nely(i)-1)',ones(nelx(i),1));
        c = repmat(b,nelz(i),1)+kron((0:nelz(i)-1)',...
            nonx*nony*ones(nelx(i)*nely(i),1));
        d = repmat((1:nelx(i)*nely(i))',nelz(i),4)+c;
        eleNode = [d,d+nonx*nony];
        if i>1
            eleNode = nodes(eleNode);
        end
        ElNodes = [ElNodes;eleNode]; %ok
    end
%==========================================================================
    % EDofMat = Degrees of freedom for each element
    EDofMat = kron(ElNodes,[3,3,3])+repmat([-2,-1,0],1,8);
%==========================================================================
    % noe = Number of elements
    % ndof = Number of degrees of freedom
    % nonodes = Number of nodes
    noe = sum(nelx.*nely.*nelz);
    nonodes = max(ElNodes,[],'all');
    ndof = 3*nonodes;
%==========================================================================
    % NodeRep = reppetition of each node in all element;
    NodeRep=groupcounts(ElNodes(:));
end
%% Boundary conditions implementation function
function [F,Dis,FixedNodes,ukdis,kdis] = ...
    BoundaryConditionsImplementation(ndof,EBC,NBC,NodeCoord,els)
%==========================================================================
    % Dis = Displacement vector
    % F = Force vector
    Dis = nan(ndof,1);
    F = zeros(ndof,1);
%==========================================================================
    for s=1:2
        if s==1
            BC = EBC;
            ND = Dis;
        else
            BC = NBC;
            ND = F;
        end
        for i=1:size(BC,1)
            Cond = NodeCoord(:,1)>=BC(i,1) & ...
                   NodeCoord(:,1)<=BC(i,2) & ...
                   NodeCoord(:,2)>=BC(i,3) & ...
                   NodeCoord(:,2)<=BC(i,4) & ...
                   NodeCoord(:,3)>=BC(i,5) & ...
                   NodeCoord(:,3)<=BC(i,6);
            Cond = find(Cond);
            if ~isempty(Cond)
                ND(Cond*3-2) = BC(i,7);
                ND(Cond*3-1) = BC(i,8);
                ND(Cond*3) = BC(i,9);
            else
                distance = sum((NodeCoord-[BC(i,1),BC(i,3),BC(i,7)]).^2,2);
                Cond = find(distance<=els);
                ND(Cond*3-2) = BC(i,7)/numel(Cond);
                ND(Cond*3-1) = BC(i,8)/numel(Cond);
                ND(Cond*3) = BC(i,9)/numel(Cond);
            end
        end
        if s==1
            Dis = ND;
        else
            F = ND;
%==========================================================================
            % FixedNodes = Nodes should remain in design 
            % domain in optimization loop, for example NBC nodes
            FixedNodes = Cond;
        end
    end
%==========================================================================
    % ukdis = Unknown displacement index
    % kdis = Known displacement index
    ukdis = isnan(Dis);
    kdis = ~isnan(Dis);
end
%% Brick element stiffness calculation function
function [Ke,Vole,Bmat,Emat] = BrickElementStiffnessCalculation(ngp,els,v)

% Notes!!!
% This function is usable for structres with elements with equal size at 
% each direction
%==========================================================================
% ngp = Number of gauss points for each direction
% els = Element size at each direction
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,els,0,0,els,els,0]';
    eYcor = [0,0,els,els,0,0,els,els]';
    eZcor = [0,0,0,0,els,els,els,els]';
%==========================================================================
 % Elasticity matrix by cosidering 1 as modulus of elasticity
Emat = 1/((1+v)*(1-2*v))*...
    [1-v,v,v,0,0,0;...
    v,1-v,v,0,0,0;...
    v,v,1-v,0,0,0;
    (1-2*v)/2*[zeros(3),eye(3)]];
%==========================================================================
    % Strain-displacement matrix
    Bmat = zeros(6,24);
%==========================================================================
% Volume at gussian points
    Volg = zeros(ngp,ngp,ngp);
%==========================================================================
    % Vole = Volume for each element
    Vole = 0;
%==========================================================================
    % Ke = stiffness matrix for element
    Ke = zeros(24);
%==========================================================================
% gauss points and their wheights
    [gp,wgp]=makegaussianpoint(ngp);
%==========================================================================
   % Loop over gauss points for Stiffness calcualtion
     for i=1:ngp
        t = gp(i);
        for j=1:ngp
            s = gp(j);
            for k=1:ngp
                r = gp(k);
%==========================================================================
% 8-node shape functions for brick element
%                 N=0.125*[(1-r)*(1-s)*(1-t),(1+r)*(1-s)*(1-t),...
%                          (1+r)*(1+s)*(1-t),(1-r)*(1+s)*(1-t),...
%                         (1-r)*(1-s)*(1+t),(1+r)*(1-s)*(1+t),...
%                         (1+r)*(1+s)*(1+t),(1-r)*(1+s)*(1+t)];
%==========================================================================
                % Calculation of shape function derivative wrt local
                % coordiantes
                N_r = 0.125*[-(1-s)*(1-t),(1-s)*(1-t),...
                             (1+s)*(1-t),-(1+s)*(1-t),...
                             -(1-s)*(1+t),(1-s)*(1+t),...
                            (1+s)*(1+t),-(1+s)*(1+t)];
    
                N_s = 0.125*[-(1-r)*(1-t),-(1+r)*(1-t),...
                               (1+r)*(1-t),(1-r)*(1-t),...
                             -(1-r)*(1+t),-(1+r)*(1+t),...
                              (1+r)*(1+t),(1-r)*(1+t)];
    
                N_t=0.125*[-(1-r)*(1-s),-(1+r)*(1-s),...
                           -(1+r)*(1+s),-(1-r)*(1+s),...
                            (1-r)*(1-s),(1+r)*(1-s),...
                            (1+r)*(1+s),(1-r)*(1+s)];
%==========================================================================
                % Calcualtion derivitave of global coordiantes 
                % wrt local coordiantes
                X_r = N_r*eXcor; Y_r = N_r*eYcor; Z_r = N_r*eZcor;
                X_s = N_s*eXcor; Y_s = N_s*eYcor; Z_s = N_s*eZcor;
                X_t = N_t*eXcor; Y_t = N_t*eYcor; Z_t = N_t*eZcor;
%==========================================================================
                % Jacobian matrix calculation
                J = [X_r,Y_r,Z_r;X_s,Y_s,Z_s;X_t,Y_t,Z_t];
                N_X_Y_Z = J\[N_r;N_s;N_t];
%==========================================================================
                % For Normal strain at x direction
                    Bmat(1,1:3:end)=N_X_Y_Z(1,:);

                % For Normal strain at y direction
                    Bmat(2,2:3:end)=N_X_Y_Z(2,:);

                % For Normal strain at z direction
                    Bmat(3,3:3:end)=N_X_Y_Z(3,:); 

                % For Shear strain at x and y direction
                    Bmat(4,1:3:end)=N_X_Y_Z(2,:);
                    Bmat(4,2:3:end)=N_X_Y_Z(1,:);

                % For Shear strain at x and z direction
                    Bmat(5,1:3:end)=N_X_Y_Z(3,:);
                    Bmat(5,3:3:end)=N_X_Y_Z(1,:);

                % For Shear strain at y and z direction
                    Bmat(6,2:3:end)=N_X_Y_Z(3,:);
                    Bmat(6,3:3:end)=N_X_Y_Z(2,:);
%==========================================================================
                % Volume at gussian points
                    Volg(k,j,i) = det(J)*wgp(k)*wgp(j)*wgp(i);

                % Volume for each element = sum(Arg)
                    Vole = Vole+Volg(k,j,i);
%==========================================================================
                % Stiffness for each element = sum(Stiffness at each
                % gussian point)
                    Ke = Ke+Bmat'*Emat*Bmat*Volg(k,j,i);
%==========================================================================
            end
        end
    end
end
%% T1 and T2 Generator function
function [T1e,T2e] = RD_T1T2Generator(ngp,els)

% Notes!!!
% This function is usable for structres with elements with equal size at 
% each direction
%==========================================================================
% ngp = Number of gauss points for each direction
% els = Element size at each direction
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,els,0,0,els,els,0]';
    eYcor = [0,0,els,els,0,0,els,els]';
    eZcor = [0,0,0,0,els,els,els,els]';
    % Volume at gussian points
    Volg = zeros(ngp,ngp,ngp);
%==========================================================================
    % In RDLS 
    % T1 = integral(N'*N)
    % T2 = integral(grad(N')*grad(N))
    % where N is shape function vector with the size of 1X8 and grad(N) and
    % grad(N) is gradient of shape function vector wrt x,y and z with the
    % size of 3X8
    T1e = zeros(8);
    T2e = zeros(8);
%==========================================================================
% gauss points and their wheights
    [gp,wgp]=makegaussianpoint(ngp);
%==========================================================================
   % Loop over gauss points for Stiffness calcualtion
     for i=1:ngp
        t = gp(i);
        for j=1:ngp
            s = gp(j);
            for k=1:ngp
                r = gp(k);
%==========================================================================
% 8-node shape functions for brick element
                N=0.125*[(1-r)*(1-s)*(1-t),(1+r)*(1-s)*(1-t),...
                         (1+r)*(1+s)*(1-t),(1-r)*(1+s)*(1-t),...
                        (1-r)*(1-s)*(1+t),(1+r)*(1-s)*(1+t),...
                        (1+r)*(1+s)*(1+t),(1-r)*(1+s)*(1+t)];
%==========================================================================
                % Calculation of shape function derivative wrt local
                % coordiantes
                N_r = 0.125*[-(1-s)*(1-t),(1-s)*(1-t),...
                             (1+s)*(1-t),-(1+s)*(1-t),...
                             -(1-s)*(1+t),(1-s)*(1+t),...
                            (1+s)*(1+t),-(1+s)*(1+t)];
    
                N_s = 0.125*[-(1-r)*(1-t),-(1+r)*(1-t),...
                               (1+r)*(1-t),(1-r)*(1-t),...
                             -(1-r)*(1+t),-(1+r)*(1+t),...
                              (1+r)*(1+t),(1-r)*(1+t)];
    
                N_t=0.125*[-(1-r)*(1-s),-(1+r)*(1-s),...
                           -(1+r)*(1+s),-(1-r)*(1+s),...
                            (1-r)*(1-s),(1+r)*(1-s),...
                            (1+r)*(1+s),(1-r)*(1+s)];
%==========================================================================
                % Calcualtion derivitave of global coordiantes 
                % wrt local coordiantes
                X_r = N_r*eXcor; Y_r = N_r*eYcor; Z_r = N_r*eZcor;
                X_s = N_s*eXcor; Y_s = N_s*eYcor; Z_s = N_s*eZcor;
                X_t = N_t*eXcor; Y_t = N_t*eYcor; Z_t = N_t*eZcor;
%==========================================================================
                % Jacobian matrix calculation
                J = [X_r,Y_r,Z_r;X_s,Y_s,Z_s;X_t,Y_t,Z_t];
                N_X_Y_Z = J\[N_r;N_s;N_t];
%==========================================================================
                % Volume at gussian points
                    Volg(k,j,i) = det(J)*wgp(k)*wgp(j)*wgp(i);
%==========================================================================
                % T1 and T2 for each element = sum(Stiffness at each
                % gussian point)
                    T1e = T1e+N'*N*Volg(k,j,i);
                    T2e = T1e+N_X_Y_Z'*N_X_Y_Z*Volg(k,j,i);
%==========================================================================
            end
        end
    end
end
