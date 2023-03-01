clear all
close all
format compact

disp(' ')
syms x y z w t
global uex uu


% THIS SCRIPT WORKS FOR THE 2D CASE ONLY (1D space + time OR 2D space)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------- stencil ------------------
% Choice of the stencil to discretize first and second order derivatives
% and respective accuracy order
stencil_First_ord=-[-1 0 1]/2; ord2=2;  % ( -u_{j-1}+u_{j+1} )/2
stencil_Laplacian=[1 -2 1]; ord1=2;     % u_{j-1}-2u_{j}+u_{j+1}
% stencil_First_ord=-[-1 1 0]/2; ord2=1;


MAXIT=20; % maximum amount of iterations 
N_c=5;    % number of coarse mesh points in each dimension (N_c in paper)
N_f=100;  % number of fine mesh points in each dimension (N_f in paper)
corr=3; % correction: 1=siolution dense in x is more accurate, 2= solution dense in y is more accurate, 3=pick randomly

% settings of the PDE
d1=1+x.^2;    d2=2+x.*y;   %diffusion coefficients (alpha_1,alpha_2 in the paper)
nu1=2-x;      nu2=1+y;     %velocity (beta_1,beta_2 in the paper)

d1=1+0*x; d2=1+0*x; n1=1+0*x;nu2=1+0*x;

TIME=0; % 0=no time dependency, 1=time-dependency
Tdiscr=1; % time discretization method: 1= Crank Nicolson, 2= Implicit Euler
DIM=2-TIME;  % space dimension


if TIME+DIM~=2
    error('This script only works for 2D cases (considering time as additional dimension)')
end

% True solution:
if TIME==0
    uu=sin(x*pi+y*pi); % solution to the PDE (use dot product if needed)
else
    uu=sin(x*pi+t*pi); % solution to the PDE (use dot product if needed)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF THE ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N{1} = [N_f,N_c*ones(1,DIM-1),N_c*ones(1,TIME)];
Coarse=[N_c*ones(1,DIM)-1,N_c*ones(1,TIME)];
Fine=[N_f*ones(1,DIM)-1,N_f*ones(1,TIME)];


space_ord=min(ord1,ord2); %accuracy order



a1=0; a2=0; str=[];
for i=1:(DIM+TIME+1)
    if i<=DIM
        if i==1
            axis='x';
        elseif i==2
            axis='y';
            N{i}=circshift(N{i-1},1);
        elseif i==3
            axis='z';
            N{i}=circshift(N{i-1},1);
        elseif i==4
            axis='w';
            N{i}=circshift(N{i-1},1);
        end
    end
    if i==(DIM+TIME) && TIME==1
        axis='t';
        N{i}=circshift(N{i-1},1);
    end
    if i==(DIM+TIME+1) %needed for richardson extrapolation
        axis='COARSE';
        N{i}=N_c*ones(1,DIM+TIME);
    end
        
    Nx=N{i}(1)-1;
    if DIM>1
        Ny=N{i}(2)-1;
    end
    if DIM>2
        Nz=N{i}(3)-1;
    end
    if DIM>3
        Nw=N{i}(4)-1;
    end
    if TIME==1
        Nt=N{i}(end);
    end



    u=char(uu);
    if TIME==1
        if DIM==1
            d1=d1+0*x; nu1=nu1+0*x;
            var{i}={'x',Nx;'t',Nt};
            f=char(diff(uu,1,'t')-(d1.*diff(uu,2,'x')-nu1.*diff(uu,1,'x')));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                [],f,f};
        end
    else
        if DIM==2
            d1=d1+0*x; d2=d2+0*x; nu1=nu1+0*x;nu2=nu2+0*x;
            var{i}={'x',Nx;'y',Ny};
            f=char(-d1.*diff(uu,2,'x')-d2.*diff(uu,2,'y')+nu1.*diff(uu,1,'x')+nu2.*diff(uu,1,'y'));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'y',char(d2),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                'y',char(nu2),@(h)stencil_First_ord/h;
                [],f,f};    
        end
    end

    if TIME==1
        if Tdiscr==1
            LH={'CN'}; ord_time=2;
        elseif Tdiscr==2
            LH={'IE'}; ord_time=1;
        else
            LH={'EE'}; ord_time=1;
        end
        if i==1
            fprintf('Time discretization = %s   -   ',LH{1})
        end
    else
        ord_time=0;
        LH={};
        if i==1
            fprintf('No Time   -   ')
        end
    end
    if i==1
        fprintf('(N_c = %.0f,  N_f = %.0f )\r',N_c,N_f)
    end
    tic
    [Ax,iAx{i},b{i},uex{i},grid{i},A{i}] = LinearSystemGenerator(var{i},LH,RH,u);
    a1=a1+toc;
    tic 
    U{i}=reshape(iAx{i}(b{i}),grid{i}.GP);
    
    err(i)=norm(uex{i}(:)-U{i}(:))/norm(uex{i}(:));
    if i<3
        str=[str,'Error mesh dense in ',axis,': %0.2e   -   '];
    else
        str=[str,'Error ',axis,' mesh: %0.2e   -   '];
    end
    h{i}=grid{i}.h;
    dt{i}=grid{i}.dt;
    a2=a2+toc;
end

%generate fine linear system
tic
var2=var{1};
for i=1:size(var{1},1)
    var2{i,2}=var{i}{i,2};
%     Fine(i)=var{i}{i,2};
end
[Ax_fine,iAx_fine,b_fine,uex_fine,grid_fine,A_fine] = LinearSystemGenerator(var2,LH,RH,u);
a1=a1+toc;


e_min=norm(iAx_fine(b_fine(:))-uex_fine(:))/norm(uex_fine(:));
    
fprintf('Time needed to generate the linear systems = %0.1fs, for the solution = %.1fs\r',a1,a2)
fprintf('Minimum possible error (error given by the fine linear system) : %.2e\r',e_min)
fprintf([str,'\r'],err)



% Keep only the common coarse points (coarse structure)
for i=1:length(U)
    P{i}=Projector(grid{i}.GP,Coarse,LH);
    P_bc{i}=ProjectorBC(Fine,grid{i}.GP,LH);
    Pf2a{i}=Projector(Fine,grid{i}.GP,LH);
    Uc{i}=P{i}*U{i}(:);
end

% Richardson Extrapolation
ORD=[space_ord ord_time];
Uextrap=RichardsonExtrap(Uc,h,dt,LH,ORD);
Uextrap=reshape(Uextrap,Coarse);
err_extrap=norm(uex{end}(:)-Uextrap(:))/norm(uex{end}(:));
res_extrap=norm(b{end}(:)-A{end}*Uextrap(:))/norm(b{end}(:));

fprintf('Error COARSE after extrapolation = %0.2e,  Residual = %0.1e     \r',err_extrap,res_extrap)



% Correction
tic
[Ucorr,U]=Correction(U,P,Pf2a,Uextrap,grid);
t_corr=toc;

Ucorr=reshape(Ucorr,Fine);
% for i=1:length(U)
%     e(i)=norm(U{i}(:)-uex{i}(:))/norm(uex{i}(:));
% end
% fprintf(['Error after correction : \r',str,' required time = %.1fs\r'],e,t_corr) %norm(Uup(:)-Uf(:))/norm(Uf(:))



% FILLING HOLES (UPSAMPLING)
tic
Ucorr=FillingHoles(Ucorr,Fine,Coarse,LH,A_fine,b_fine,TIME);

Ucorr=reshape(Ucorr,Fine);
e_up=norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:));
r_up=norm(b_fine(:)-A_fine*Ucorr(:))/norm(b_fine(:));
t_up=toc;
fprintf(['Error after UPsampling           = %0.2e,  Residual = %0.1e     (solution time = %.1fs)\r'],e_up,r_up,t_up) %norm(Uup(:)-Uf(:))/norm(Uf(:))





% ITERATE THE ALGORITHM
clear e
fig=0;




rich=0; % cannot use richardson extrapolation for the error equation
ORD=[space_ord ord_time];
for k=1:MAXIT
    res=b_fine(:)-A_fine*Ucorr(:);
    errr(k)=norm(uex_fine(:)-Ucorr(:))/norm(uex_fine(:));
    res_norm(k)=norm(res)/norm(b_fine);
    fprintf(['IT=%.0f - ERROR : %0.3e      RESIDUAL : %0.3e\r'],k,norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:)),norm(res)/norm(b_fine(:))) %norm(Uup(:)-Uf(:))/norm(Uf(:))

    
    b_fine=b_fine(:); res=res(:);
        
    e={};
%     e_ex=A_fine\b_fine;
%     PP=Pf2a;
    PP=P_bc;
    for proj=1:length(P)
%         ee_ex{proj}=PP{proj}*e_ex;
        r{proj}=PP{proj}*res;
        e{proj}=A{proj}\r{proj};
%         norm(e{proj})/norm(ee_ex{proj})
        %normalization:
        if proj==1
            nE=norm(e{proj});
        else
            e{proj}=e{proj}/norm(e{proj})*nE;
        end
    end
    

    % Richardson Extrapolation
    rich=0;
    if rich==1
        for proj=1:length(P)
            ec{proj}=P{proj}*e{proj};   % P{proj}
        end
        Uextrap=RichardsonExtrap(ec,h,dt,LH,ORD);
        Uextrap=reshape(Uextrap,Coarse);
        fprintf(['Accuracy of projected error (before extrapolation) : \r',str,'\r'],e2) %norm(Uup(:)-Uf(:))/norm(Uf(:))
        
        err_extrap=norm(Pf2a{end}*e_ex(:)-Uextrap(:))/norm(Pf2a{end}*e_ex(:));
        fprintf('Accuracy of COARSE error after extrapolation = %0.2e\r',err_extrap)
    else
        if corr==3
            corrr=1+round(rand);
            Uextrap=P{corrr}*e{corrr};
        else
            Uextrap=P{corr}*e{corr};
        end
    end

    
    % Correction
    [Ecorr,~]=Correction(e,P,Pf2a,Uextrap,grid);

    
    % Filling holes
    Ecorr=FillingHoles(Ecorr,Fine,Coarse,LH,A_fine,zeros(size(b_fine)),TIME);
    
    %Apply correction
    f=@(c)norm(b_fine(:)-A_fine*reshape(Ucorr(:)+c*Ecorr(:),[],1));
    OPTIONS.TolX=10^-10;
    cc2=fminbnd(f,0,10,OPTIONS);
    Ucorr=Ucorr(:)+cc2*Ecorr;
    e{end+1}=cc2*Ecorr;
    
    
end



fprintf(['IT=%.0f - ERROR : %0.3e      RESIDUAL : %0.3e\r'],k,norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:)),norm(res)/norm(b_fine(:))) %norm(Uup(:)-Uf(:))/norm(Uf(:))

return


% z=abs(reshape(res,N_f-1,N_f-1));
% x=grid_fine.X{1}(2:end-1,2:end-1);
% y=grid_fine.X{2}(2:end-1,2:end-1);
% scatter(x(:),y(:),7,z(:)+eps,'filled')
% colorbar()


%time
z=abs(reshape(res,N_f-1,N_f));
x=grid_fine.X{1}(2:end-1,1:end-1);
y=grid_fine.T(2:end-1,1:end-1);
surf(x,y,z)

%space
z=abs(reshape(res,N_f-1,N_f-1));
x=grid_fine.X{1}(2:end-1,2:end-1);
y=grid_fine.X{2}(2:end-1,2:end-1);
surf(x,y,z)














function Ucorr=FillingHoles(Ucorr,Fine,Coarse,LH,A_fine,b_fine,TIME)
P_f2holes=Projector_Fine2Holes(Fine,Coarse,LH);
A2=P_f2holes*(A_fine*P_f2holes');
f2=P_f2holes*b_fine(:);

% blocks=A2\(f2-P_f2holes*A_fine*Ucorr(:));
for i=1:length(Fine)
    if i==length(Fine) && TIME==1
        hole(i)=Fine(i)/Coarse(i)-1;
        space(i)=hole(i)*(Coarse(i));
    else
        hole(i)=(Fine(i)+1)/(Coarse(i)+1)-1;
        space(i)=hole(i)*(Coarse(i)+1);
    end
end

b=f2-P_f2holes*(A_fine*Ucorr(:));
blocks=zeros(size(b));

if TIME==1
    B1=A2(1:prod(hole),1:prod(hole));
    for it=1:2
        blocks=blocks+reshape(B1\reshape(b-A2*blocks,prod(hole),[]),[],1);
    end
    [blocks,FLAG,RELRES,ITER,RESVEC]=gmres(A2,b,100,10^-14,[],@(r)reshape(B1\reshape(r,prod(hole),[]),[],1),[],blocks);
else
    v=hole(1)*prod(space(2:end));
    B1=A2(1:v,1:v);
    for it=1:2
        blocks=blocks+reshape(B1\reshape(b-A2*blocks,v,[]),[],1);
    end
    [blocks,FLAG,RELRES,ITER,RESVEC]=gmres(A2,f2-P_f2holes*A_fine*Ucorr(:),100,10^-14,[],@(r)reshape(B1\reshape(r,v,[]),[],1),[],blocks);
end
Ucorr=Ucorr(:)+P_f2holes'*blocks;
end





function P=Projector_Fine2Holes(Fine,Coarse,LH)
% amount of grid points in vector "Coarse" and "Fine" are ordered according
% to the ordering in "LinearSystemGenerator"
    P=1;
    for i=1:length(Fine)
        if Fine(i)==Coarse(i)
            P=kron(speye(Coarse(i)),P);
        else
            if (i==length(Fine)) && ~isempty(LH) %equation depends on time (keep last time step)
                PP=speye(Fine(i));
                s=Fine(i)/Coarse(i);
                PP(s:s:Fine(i),:)=[];
            else
                PP=speye(Fine(i));
                s=(Fine(i)+1)/(Coarse(i)+1);
                PP(s:s:Fine(i)-1,:)=[];
            end
            P=kron(PP,P);
        end
    end
end






function [Ucorr,U]=Correction(U,P,Pf2a,Uextrap,grid)
for i=1:length(U)
    Err{i}=reshape(Uextrap,grid{end}.GP)-reshape(P{i}*U{i}(:),grid{end}.GP);
    Uu{i}=reshape(P{i}'*Uextrap(:),grid{i}.GP);
end

ord=5;
for kk=1:length(grid{end}.xx)
%     xx{kk}=grid{end}.xx{kk};
    ORD{kk}=ord;
end
if isempty(grid{end}.tt)==0
%     xx{kk+1}=grid{end}.tt;
    ORD{kk+1}=ord;
end
xx=[grid{end}.xx,grid{end}.tt];
for k=1:size(U,2) % dimensions
%     e=Err{k};%reshape(P{k}*Err{k}(:),grid{end}.GP);
%     p=spapi(ord,[0;grid{end}.xx',1],[0;e(i,:)']); %spline interp
    e=padarray(Err{k},ones(1,length(grid{end}.GP)),0);
    if isempty(grid{end}.tt)==0
        if length(size(U{1}))==2
            e=e(:,1:end-1);
        elseif length(size(U{1}))==3
            e=e(:,:,1:end-1);
        elseif length(size(U{1}))==4
            e=e(:,:,:,1:end-1);
        elseif length(size(U{1}))==5
            e=e(:,:,:,:,1:end-1);
        end
    end
    p=spapi(ORD,xx,e); %spline interp
    U{k} = reshape(U{k},grid{k}.GP)+fnval(p,[grid{k}.xx_intern,grid{k}.tt_intern]);
    if k==1
    	Ucorr=Pf2a{k}'*U{k}(:); % Projector fine to anisotropic mesh
    else
        Ucorr=Ucorr+Pf2a{k}'*U{k}(:)-Pf2a{k}'*(P{k}'*(P{k}*U{k}(:)));
    end
end
end






% profile viewer
function Uextrap=RichardsonExtrap(U,h,dt,LH,ORD)
% Removes lowest order terms
for i=1:length(U)-1
    for j=1:length(h{1})
        LS(j,i)=h{i}(j)^ORD(1);% generating linear system
    end
    if ~isempty(LH) %equation depends on time
        if strcmp(LH{1},'CN') % order 2 accurate
            LS(j+1,i)=dt{i}^ORD(2);
        else %IE and EE are of order 1
            LS(j+1,i)=dt{i}^ORD(2);
        end
    end
    if i<=length(h{end})
        b(i,1)=h{end}(i)^ORD(1);
    end
end
if ~isempty(LH) %equation depends on time
    if strcmp(LH{1},'CN') % order 2 accurate
        b(end+1,1)=dt{end}^ORD(2);
    else %IE and EE are of order 1
        b(end+1,1)=dt{end}^ORD(2);
    end
else
%     b(i,1)=h{1}(end)^ORD(1);
end

sol=LS\(-b);
sol=[sol;1];
Uextrap=0;
for i=1:length(U)
    Uextrap=Uextrap+sol(i)*U{i}(:);
end
Uextrap=Uextrap/sum(sol);
end



function P=Projector(Fine,Coarse,LH)
% amount of grid points in vector "Coarse" and "Fine" are ordered according
% to the ordering in "LinearSystemGenerator"
P=1;
for i=1:length(Fine)
    if Fine(i)==Coarse(i)
        P=kron(speye(Coarse(i)),P);
    else
        if (i==length(Fine)) && ~isempty(LH) %equation depends on time (keep last time step)
            PP=speye(Fine(i));
            s=Fine(i)/Coarse(i);
            PP=PP(s:s:Fine(end),:);
        else
            PP=speye(Fine(i));
            s=(Fine(i)+1)/(Coarse(i)+1);
            PP=PP(s:s:Fine(i)-1,:);
        end
        if size(PP,1)~=Coarse(i)
            error('something is wrong')
        end
        P=kron(PP,P);
    end
end
end








function P=ProjectorBC(Fine,Coarse,LH)
% amount of grid points in vector "Coarse" and "Fine" are ordered according
% to the ordering in "LinearSystemGenerator"
P=1;
for i=1:length(Fine)
    if Fine(i)==Coarse(i)
        P=kron(speye(Coarse(i)),P);
    else
        if (i==length(Fine)) && ~isempty(LH) %equation depends on time (keep last time step)
            if strcmp(LH{1},'CN')
                PP=speye(Fine(i));
                s=Fine(i)/Coarse(i);
                PP=PP(s:s:Fine(end),:);
                for j=1:size(PP,1)-1
                    PP(j+1,s*j)=1;
                end
                PP(1,1)=1;
                PP=PP/2;
%                 
            elseif strcmp(LH{1},'IE')
                PP=speye(Fine(i));
                s=Fine(i)/Coarse(i);
                PP=PP(s:s:Fine(end),:);
            elseif strcmp(LH{1},'EE')
                error('not implemented yet')
            end
        else
            PP=speye(Fine(i));
            s=(Fine(i)+1)/(Coarse(i)+1);
            PP=PP(s:s:Fine(i)-1,:);
        end
        if size(PP,1)~=Coarse(i)
            error('something is wrong')
        end
        P=kron(PP,P);
    end
end
end




