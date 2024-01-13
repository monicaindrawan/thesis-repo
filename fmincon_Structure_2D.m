function history = fmincon_top88_2D
tic
%% SET UP SHARED VARIABLES WITH OUTFUN
history.x = [];
history.fval = [];
%% MATERIAL PROPERTIES
nely = 30;
nelx = 90;
volfrac = 0.5;
penal = 3.0;
rmin = 2.5;
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% SOLUTION INITIALIZATION
% VOLUME FRACTION
x0(1:nelx,1:nely) = volfrac;
xm = x0(:);
eq = ones(1,nelx*nely);
totalvolfrac = nely*nelx*volfrac;
% UPPER AND LOWER BOUNDS
nvar=nelx*nely;
    lb(1:nvar) = 0.0;
    ub(1:nvar) = 1.0;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS
% MBB PROBLEM
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% CANTILEVER PROBLEM
F = sparse(2*(nely+1)*(nelx+1),1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% BRIDGE PROBLEM
F = sparse((nely+1)*(nelx)+2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [(2*nely)+1,2*(nely+1),2*(nely+1)*(nelx+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% CALLBACK OBJECTIVE FUNCTION AND GRADIENT OBJECTIVE
objgrad = @(x)funcfdf2D(nelx,nely,penal,E0,Emin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,x);
hessfunc = @(x,lamda)hessinterior(nely,nelx,penal,E0,Emin,edofMat,KE,iK,jK,freedofs,U,F,nvar,x);
options = optimoptions('fmincon','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','cg','TolFun',1e-6,'TolCon',1e-8,'MaxIteration',1000,...
    'MaxFunctionEvaluations', 10000,'HessianFcn',hessfunc)
% Run FMINCON.
delete FMINCON_2D_Data.txt
[x,fval] = fmincon(objgrad,xm,[],[],eq,totalvolfrac,lb,ub,[],options);
FMINCONtime = toc;
save('FMINCON_2D_Time.mat','FMINCONtime');
% HISTORICAL DATA COMPLIANCE
function stop = outfun(x,optimValues,state)
     stop = false;
     switch state
         case 'iter'
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x x];
     end
end
%% PLOT DENSITIES
iter = size(history.x);
 for i = 1:iter(2)
     x2 = vec2mat(history.x(:,i)',nelx);
     xmat = 1-x2;
     colormap(gray); imagesc(1-x2); caxis([0 1]); axis equal; axis off; drawnow;
 end
 save('FMINCON_2D_Topology.mat','xmat');
end
% ----------------------------------------------------------------------
function [c,dcvect] = funcfdf2D(nelx,nely,penal,E0,Emin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,x)
warning('off','all');
xPhys = vec2mat(x,nelx);
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
%% FILTERING/MODIFICATION OF SENSITIVITIES
dc(:) = H*(xPhys(:).*dc(:))./Hs./max(1e-3,xPhys(:));
dctrans = dc';
dcvect = dctrans(:);
%% SAVING COMPLIANCE DATA AND CPU TIME
FEtime = toc;
diary FMINCON_2D_Data.txt
fprintf("%f",c)
fprintf("\n")
diary off 
end
% ----------------------------------------------------------------------
function Hess = hessinterior (nely,nelx,penal,E0,Emin,edofMat,KE,iK,jK,freedofs,U,F,nvar,x)
xPhys = vec2mat(x,nelx);
% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
% SECOND DERIVATIVE
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
ddc = 2*(penal*(E0-Emin)*xPhys.^(penal-1)).^2.*(1./(Emin+xPhys.^penal*(E0-Emin))).*ce;

Hvect = zeros(nvar,1); a = zeros(nvar,1); m = 1;
for i = 1:nely
    for j = 1:nelx
            a(m) = m;
            Hvect(m) = ddc(i,j);
            m = m + 1;
        end
end
Hess = sparse(a,a,Hvect,nvar,nvar);
end