function history = fmincon_top88_3D
tic
%% SET UP SHARED VARIABLES WITH OUTFUN
history.x = [];
history.fval = [];
%% MATERIAL PROPERTIES
nelx = 40;
nely = 40;
nelz = 5;
volfrac = 0.3;
rmin = 1.4;
penal = 3.0;
k0 = 1;
kmin = 1e-3;
%% SOLUTION INITIALIZATION
% VOLUME FRACTION
x0(1:nely,1:nelx,1:nelz) = volfrac;
xm = x0(:);
eq = ones(1,nele);
totalvolfrac = nele*volfrac;
% UPPER AND LOWER BOUNDS 
nele=nelx*nely*nelz;
    lb(1:nele) = 0.0;
    ub(1:nele) = 1;
%% DEFINE HEAT GENERATION AND BOUNDARY CONDITION
il = nelx/2-nelx/20:nelx/2+nelx/20; jl = nely; kl = 0:nelz;                          % Coordinates
fixedxy = il*(nely+1)+(nely+1-jl);                                                   % Coordinates
fixednid = repmat(fixedxy',size(kl))+repmat(kl*(nelx+1)*(nely+1),size(fixedxy,2),1); % Node IDs
fixeddof = reshape(fixednid,[],1);                                                   % DOFs
%% PREPARE FINITE ELEMENT ANALYSIS
ndof = (nelx+1)*(nely+1)*(nelz+1);
F = sparse(1:ndof,1,-0.01,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(k0);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = nodeids(:)+1;
edofMat = repmat(edofVec,1,8)+ ...
    repmat([0 nely + [1 0] -1 ...
    (nely+1)*(nelx+1)+[0 nely + [1 0] -1]],nele,1);
iK = reshape(kron(edofMat,ones(8,1))',8*8*nele,1);
jK = reshape(kron(edofMat,ones(1,8))',8*8*nele,1);
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
objgrad = @(x)funcfdf3D(nely,nelx,nelz,penal,k0,kmin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,nele,x);
hessfunc = @(x,lamda)hessinterior(nely,nelx,nelz,penal,k0,kmin,edofMat,KE,iK,jK,freedofs,U,F,nele,x);
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
for m = 1:iter(2)-1
Rho = zeros(nely,nelx,nelz);
counter = 1;
for i = 1:nely
    for j = 1:nelx
        for k = 1:nelz
            Rho(i,j,k) = x(counter)
            counter = counter + 1;
        end
    end
end
%% PLOT DENSITIES
clf; display_3D(Rho)
end
save('FMINCON_3D_Topology.mat','Rho')
end
% === GENERATE THERMAL CONDUCTIVITY MATRIX ===
function [KE] = lk_H8(k)
A1 = 4*eye(2); A2 = -eye(2);
A3 = fliplr(A2); A4 = -ones(2);
KE1 = [A1 A2; A2 A1];
KE2 = [A3 A4; A4 A3];
KE = 1/12*k*[KE1 KE2; KE2 KE1];
end
% ----------------------------------------------------------------------
function [c,dcvect] = funcfdf3D(nely,nelx,nelz,penal,k0,kmin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,nele,x)
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
sK = reshape(KE(:)*(kmin+(1-kmin)*xPhys(:)'.^penal),8*8*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
c = sum(sum(sum((kmin+(1-kmin)*xPhys.^penal).*ce)));
dc = -penal*(1-kmin)*xPhys.^(penal-1).*ce;
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
function Hess = hessinterior (nely,nelx,nelz,penal,k0,kmin,edofMat,KE,iK,jK,freedofs,U,F,nele,x)
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
sK = reshape(KE(:)*(kmin+(1-kmin)*xPhys(:)'.^penal),8*8*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
%% SECOND DERIVATIVE
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
ddc = 2*(penal*(1-kmin)*xPhys.^(penal-1)).^2.*(1/(kmin+(1-kmin)*xPhys.^penal)).*ce;
%% HESSIAN MATRIX
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
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
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