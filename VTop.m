function [x,vem] = VTop(vem,opt)
%   VTOP  a 2D proportional topology optimization (PTO) method using
%   virtual elements. This code was written basing on the structure of
%   Polytop by Talischi et al.(https://doi.org/10.1007/s00158-011-0696-x)
%   Input
%   vem:    parameters setup of the domain using virtual elements
%   opt:    parameters setup of the PTO
%   Output
%   x:      material distributions
%   ----------------------------------------------------------------
%   Please cite following paper when you use this code:
%   An enhanced proportional topology optimization with virtual
%   elements: Formulation and numerical implementation
%	----------------------------------------------------------------
%   Copyright (c) 2023
%   Authors:
%   Minh Tran.      Email: minh.tt@vgu.edu.vn
%   Minh Nguyen.    Email: nguyenngocminh6@duytan.edu.vn
%	Ho Chi Minh city, Vietnam 26th-Feb-2023
%	----------------------------------------------------------------
Iter=0; TOL=opt.Tol*(opt.xMax-opt.xMin); Change=2*TOL; x=opt.xIni; E = opt.MatIntFnc(x);
[FigHandle,FigData] = InitialPlot(vem,x);
vem = GlobalIndex(vem);
while (Change >= TOL && Iter < opt.MaxIter)
    Iter        = Iter + 1;
    [U,vem]     = VEAnalysis(vem,E);
    [ce,c]      = ObjectiveFnc(vem,U,E);
    [x,Change]  = UpdateScheme(ce,x,opt);
    E           = opt.MatIntFnc(x);
    % PRINT RESULTS
    fprintf('It.:%5i ObjFunc.:%11.4f Vol.:%7.3f Change.:%7.3f\n',...
        Iter, c, sum(x.*opt.VolElem)/sum(opt.VolElem), Change);
    set(FigHandle,'FaceColor','flat','CData',1-x(FigData)); drawnow
end
end
%% Objective function
function [ce,c]=ObjectiveFnc(vem,U,E)
temp = cumsum(U(vem.i).*vem.k.*E(vem.e).*U(vem.j));
temp = temp(cumsum(vem.ElemNDof.^2));
ce = [temp(1);temp(2:end)-temp(1:end-1)]; % element compliance
c = sum(ce); % total compliance: obj. func
end
%% Update material scheme
function [xn1,Change] = UpdateScheme(ce,xn,opt)
xTarget = sum(opt.VolElem)*opt.VolFrac;
xRemaining = xTarget;
xNew = zeros(size(xn));
C = opt.P*(ce/sum(ce));
while (xRemaining > 1e-4)
    xNew = xNew + xRemaining.*C;
    xNew = max(min(xNew,opt.xMax),opt.xMin);
    xRemaining = xTarget - sum(xNew.*opt.VolElem);
end
xn1 = opt.Alpha*xn + (1-opt.Alpha)*xNew;
Change = norm(xn1-xn)/norm(xn); % relaxed change
end
%% Finite Element Analysis
function [U,vem] = VEAnalysis(vem,E)
K = sparse(vem.i,vem.j,E(vem.e).*vem.k);
K = (K+K')/2;
% Force
NLoad = size(vem.Load,1);
vem.F = zeros(2*vem.NNode,1);  %external load vector
vem.F(2*vem.Load(1:NLoad,1)-1) = vem.Load(1:NLoad,2);% concentrated Fx
vem.F(2*vem.Load(1:NLoad,1))   = vem.Load(1:NLoad,3);% concentrated Fy
% Boundary conditions
bdx = unique(vem.Supp(vem.Supp(:,2) == 1, 1)); % nodeIDs where ux = 0
bdy = unique(vem.Supp(vem.Supp(:,3) == 1, 1)); % nodeIDs where uy = 0
fbx = bdx*2-1; % constraint in x
fby = bdy*2;   % constraint in y
FixedDofs = [fbx; fby];
BOUNDARYCOND = [FixedDofs, zeros(size(FixedDofs))];
% Solve
[~, U] = LinearSolve(K, vem.F, BOUNDARYCOND, 2*vem.NNode);
end
%% Globalize index (vectorize global index for global stiffness matrix K)
function vem = GlobalIndex(vem)
vem.ElemNDof = 2*cellfun(@length,vem.Element); % # of DOFs per element
vem.i = zeros(sum(vem.ElemNDof.^2),1);
vem.j=vem.i; vem.k=vem.i; vem.e=vem.i;
index = 0;
for eID = 1:vem.NElem
    Ke=LocalK(vem,eID);
    NDof = vem.ElemNDof(eID);
    eDof = reshape([2*vem.Element{eID}-1;2*vem.Element{eID}],NDof,1);
    I=repmat(eDof ,1,NDof); J=I';
    vem.i(index+1:index+NDof^2) = I(:);
    vem.j(index+1:index+NDof^2) = J(:);
    vem.k(index+1:index+NDof^2) = Ke(:);
    vem.e(index+1:index+NDof^2) = eID;
    index = index + NDof^2;
end
end
%% Local Element Stiffness
function [Ke] = LocalK(vem,eID)
D=vem.E0/(1-vem.Nu0^2)*[1 vem.Nu0 0;vem.Nu0 1 0;0 0 2*(1-vem.Nu0)];

eNode =   vem.Element{eID};
verts =   vem.Node(eNode,:);
area  =   vem.AreaElem(eID);
% VEM element matrices
[Wc,Hc,Wr,Hr]=Element_Projector(verts,area);
Pp=Hr*Wr'+Hc*Wc';
I_2N=eye(2*length(verts));
% scaling parameter
gamma=1.0; alpha=gamma*area*sum(diag(D))/sum(diag(Hc'*Hc));
Se=alpha*I_2N;
% VEM element stiffness
consistency = area*Wc*D*Wc';
stability   = (I_2N-Pp)'*Se*(I_2N-Pp);
Ke = vem.Thickness*(consistency + stability);
end
%-------------------------------------------------------PROJECTION OPERATOR
function [Wc,Hc,Wr,Hr] = Element_Projector(verts,area)
x_bar = (1/length(verts))*sum(verts);
n=normals_to_edges(verts);
n1 = [n(end,:); n(1:end-1,:)];
v1 = [verts(end,:); verts(1:end-1,:)];
len1 = sqrt(sum((verts-v1).^2,2));
len = len1([2:end,1]);

Wc=zeros(2*length(verts),3);
Hc=zeros(2*length(verts),3);
Wr=zeros(2*length(verts),3);
Hr=zeros(2*length(verts),3);

qa1 = 0.25*(len1.*n1(:,1)+len.*n(:,1))/area;
qa2 = 0.25*(len1.*n1(:,2)+len.*n(:,2))/area;
qa = [qa1 qa2];
nVerts = length(verts); zeroVec = zeros(nVerts,1); onesVec = ones(nVerts,1);
Wc(1:2:2*nVerts,:)=[2*qa(:,1),zeroVec,qa(:,2)]; % odd rows
Wc(2:2:2*nVerts,:)=[zeroVec,2*qa(:,2),qa(:,1)]; % even rows
Wr(1:2:2*nVerts,:)=[1/nVerts*onesVec,zeroVec,qa(:,2)]; % odd rows
Wr(2:2:2*nVerts,:)=[zeroVec,1/nVerts*onesVec,-qa(:,1)];% even rows
Hc(1:2:2*nVerts,:)=[verts(:,1)-x_bar(1), zeroVec, verts(:,2)-x_bar(2)];% odd rows
Hc(2:2:2*nVerts,:)=[zeroVec, verts(:,2)-x_bar(2), verts(:,1)-x_bar(1)];% even rows
Hr(1:2:2*nVerts,:)=[onesVec,zeroVec,verts(:,2)-x_bar(2)];% odd rows
Hr(2:2:2*nVerts,:)=[zeroVec,onesVec,-verts(:,1)+x_bar(1)];% even rows
end
function [res] = normals_to_edges(verts)
nVerts = length(verts);
v1=[verts; verts(1,:)];
w = v1(2:nVerts+1,:)-v1(1:nVerts,:);
n = [w(:,2) -w(:,1)];
% res= n./vecnorm(n,2,2);
n1= n(:,1)./sqrt(plus(n(:,1).^2,n(:,2).^2));
n2= n(:,2)./sqrt(plus(n(:,1).^2,n(:,2).^2));
res = [n1 n2];
end
%% Linear Solver
function [Fc,u,F] = LinearSolve(K, F, bc, tdof, pre_opt, LU_TOL, pre_iter)
%Purpose: impose boundary conditions and partition the stiffness matrix and
%the load vector
%--------------------------------------------------------------------------
%Written by: MINH NGUYEN
%14/03/2012, Bochum
%--------------------------------------------------------------------------
%Last reviewed and edited by: MINH NGUYEN
%14/03/2012
%--------------------------------------------------------------------------
%INPUT:       
    % K: stiffness matrix
    % F: load vector
    % bc: boundary condition
            % 1st column: index of constrained dof
            % 2nd column: value of constrained dof
    % tdof: total dofs
%--------------------------------------------------------------------------
%OUTPUT:
    % Fc: reaction force
    % u: unknown variables    
    % F: total external load
%--------------------------------------------------------------------------
% preconditioner?
if nargin < 5
    pre_opt = 0;
end
if nargin < 6
    LU_TOL = 1.e-6;
end
if nargin < 7
    pre_iter = 20;
end
% Sort bc
cdof = bc(:,1); % constrained dof
[cdof, order] = sort(cdof, 'ascend');
bc = bc(order,1:2);

udof = (1:tdof)';
udof(cdof(:)) = []; % unconstrained dof
% known values
uc = bc(:,2);
% Separate K and F
Kuu = K(udof(:), udof(:));
Kuc = K(udof(:), cdof(:));
Kcu = K(cdof(:), udof(:));
Kcc = K(cdof(:), cdof(:));

Fu = F(udof(:));
% Solve
% unknown values
if pre_opt == 0
    uu = Kuu\(Fu - Kuc*uc); % directly solved
else
    A = Kuu; b = Fu-Kuc*uc;
    [L2,U2] = luinc(sparse2(A),LU_TOL); % preconditioner by incomplete LU-decomposition
    uu = bicg(A,b,1e-15,pre_iter,L2,U2);
end
% reaction force
Fc = Kcu * uu + Kcc * uc;
% F
F(udof) = Fu;
F(cdof) = Fc;
% Rearrange
u = zeros(tdof,1);
u(udof) = uu;
u(cdof) = uc;
end
%% Visualization
function [handle,map] = InitialPlot(vem,x0)
Tri = zeros(length([vem.Element{:}])-2*vem.NElem,3);
map = zeros(size(Tri,1),1); index=0;
for el = 1:vem.NElem
    for enode = 1:length(vem.Element{el})-2
        map(index+1) = el;
        Tri(index+1,:) = vem.Element{el}([1,enode+1,enode+2]);
        index = index + 1;
    end
end
handle = patch('Faces',Tri,'Vertices',vem.Node,'FaceVertexCData',...
    1-x0(map),'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(gray); caxis([0 1]);
end