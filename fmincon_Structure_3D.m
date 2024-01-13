function history = fmincon_top88_3D
tic
%% SET UP SHARED VARIABLES WITH OUTFUN
history.x = [];
history.fval = [];
%% MATERIAL PROPERTIES
nelx = 60;
nely = 20;
nelz = 4;
volfrac = 0.3;
penal = 3.0;
rmin = 1.5;
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% SOLUTION INITIALIZATION
% VOLUME FRACTION
x0(1:nely,1:nelx,1:nelz) = volfrac;
xm = x0(:);
% UPPER AND LOWER BOUNDS
nele=nelx*nely*nelz;
    lb(1:nele) = 0.0;
    ub(1:nele) = 1;
eq = ones(1,nele);
totalvolfrac = nele*volfrac;
%% CANTILEVER SMALL SCALE PROBLEM
% top(60,20,4,0.3,3.0,1.5)
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% BRIDGE PROBLEM
% top(40,20,40,0.2,3.0,1.5)
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx/2, 0, nelz/2);               % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid([0 0 nelx nelx],[0 0 0 0],[0 nelz 0 nelz]); % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% CANTILEVER MEDIUM SCALE PROBLEM
% top(60,20,16,0.15,3.0,1.5)
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, nelz/2);                  % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);  % Node IDs
loaddof = 3*loadnid(:) - 1;                              % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% PREPARE FINITE ELEMENT ANALYSIS
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% CALLBACK OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
objgrad = @(x)funcfdf3D(nely,nelx,nelz,penal,E0,Emin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,nele,x);
hessfunc = @(x,lamda)hessinterior(nely,nelx,nelz,penal,E0,Emin,edofMat,KE,iK,jK,freedofs,U,F,nele,x);
% FMINCON OPERATOR SETTING
options = optimoptions('fmincon','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,...
    'TolFun',1e-6,'TolCon',1e-8,'MaxIteration',1000,'MaxFunctionEvaluations', 10000,...
    'HessianFcn',hessfunc);
% RUNNING FMINCON
delete FMINCON_3D_Data.txt
[x,fval] = fmincon(objgrad,xm,[],[],eq,totalvolfrac,lb,ub,[],options);
FMINCONtime = toc;
save('FMINCON_3D_Time.mat','FMINCONtime');
% HISTORICAL DATA COMPILATION
function stop = outfun(x,optimValues,state)
     stop = false;
     switch state
         case 'iter'          
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x x];
     end
end
%% CHANGE VECTOR X TO MATRIX
iter = size(history.x);
Rho = zeros(nely,nelx,nelz);
counter = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            Rho(i,j,k) = x(counter);
            counter = counter + 1;
        end
    end
end
%% PLOT DENSITIES
figure()
display_3D(Rho)
save('FMINCON_3D_Topology.mat','Rho')
end
% ----------------------------------------------------------------------
function [c,dcvect] = funcfdf3D(nely,nelx,nelz,penal,E0,Emin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,nele,x)
warning('off','all');
xPhys = zeros(nely,nelx,nelz);
counter = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            xPhys(i,j,k) = x(counter);
            counter = counter + 1;
        end
    end
end
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
dcvect = zeros(1,nely*nelx*nelz);
%% FILTERING AND MODIFICATION OF SENSITIVITIES
dc(:) = H*(xPhys(:).*dc(:))./Hs./max(1e-3,xPhys(:));
counter = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            dcvect(counter) = dc(i,j,k);
            counter = counter + 1;
        end
    end
end
FEtime = toc;
diary FMINCON_3D_Data.txt
fprintf("%f",c)
fprintf("\n")
diary off 
end
% ----------------------------------------------------------------------
function Hess = hessinterior (nely,nelx,nelz,penal,E0,Emin,edofMat,KE,iK,jK,freedofs,U,F,nele,x)
xPhys = zeros(nely,nelx,nelz);
m = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            xPhys(i,j,k) = x(m);
            m = m + 1;
        end
    end
end
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
%% SECOND DERIVATIVE
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
ddc = 2*(penal*(E0-Emin)*xPhys.^(penal-1)).^2.*(1/(Emin+xPhys.^penal*(E0-Emin))).*ce;
Hvect = zeros(nele,1); a = zeros(nele,1); m = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            a(m) = m;
            Hvect(m) = ddc(i,j,k);
            m = m + 1;
        end
    end
end
Hess = sparse(a,a,Hvect,nele,nele);
end
% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];
K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
% ----------------------------------------------------------------------
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.8)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]);
end