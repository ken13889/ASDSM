clear all
close all
format compact

profile on


% THIS SCRIPT WORKS FOR THE 3D CASE ONLY (2D space + time OR 3D space)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global uex POSS Fine Coarse discr nu1 nu2 nu3 d1 d2 d3 uu 
syms x y z w t
fprintf('\r\r')



TIME=1;
DIM=3-TIME; %space dimension
discr=1; % 1=Crank-Nicolson , 2=Implicit Euler
MAXIT=10; % Maximum amount if iterations



N_c=10; % Number of coarse mesh poitns
N_f=40; % Number of fine mesh points


% Settings of the PDE
d1=1+0*x;        d2=2+0*y;      d3=3+0*x;  % Diffusion
nu1=1+0*x;    nu2=2+0*y;   nu3=3+0*z;      % Velocity
if TIME==0 
    uu=sin(4*x*pi+4*y*pi+4*z*pi);  % True solution 3D space
else
    uu=sin(4*x*pi+4*y*pi+4*t*pi);  % True solution 2D space + Time
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF THE ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TIME==0
    N{1} = [N_c,N_f*ones(1,DIM-1)];
else
    N{1} = [N_c,N_f*ones(1,DIM-1+TIME)];
end
Coarse=[N_c*ones(1,DIM)-1,N_c*ones(1,TIME)];
Fine=[N_f*ones(1,DIM)-1,N_f*ones(1,TIME)];

stencil_Laplacian=[1 -2 1]; ord1=2;
stencil_First_ord=-[-1 0 1]/2; ord2=2;
% stencil_First_ord=-[-1 1 0]/2; ord2=1;
space_ord=min(ord1,ord2); %accuracy order


Nx=0; Ny=0; Nz=0; Nt=0;
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
    
    [var{i},LH,RH,u,ord_time]=LS(Nx,Ny,Nz,Nt,TIME,DIM,stencil_Laplacian,stencil_First_ord);

    if i==1
        if isempty(LH)==0
            fprintf('Time discretization = %s   -   ',LH{1})
        else
            fprintf('No Time   -   ')
        end
        fprintf('(N_c = %.0f,  N_f = %.0f )\r',N_c,N_f)
    end
    tic
    [Ax,iAx{i},b{i},uex{i},grid{i},A{i}] = LinearSystemGenerator(var{i},LH,RH,u);
    a1=a1+toc;
    tic 
    U{i}=reshape(iAx{i}(b{i}),grid{i}.GP);
    
    err(i)=norm(uex{i}(:)-U{i}(:))/norm(uex{i}(:));
    str=[str,'Error dense in ',axis,': %0.2e   -   '];
    h{i}=grid{i}.h;
    dt{i}=grid{i}.dt;
    a2=a2+toc;
end

%generate fine linear system
tic
var2=var{1};
for i=1:size(var{1},1)
    var2{i,2}=Fine(i);
%     Fine(i)=var{i}{i,2};
end
[Ax_fine,iAx_fine,b_fine,uex_fine,grid_fine,A_fine] = LinearSystemGenerator(var2,LH,RH,u);
a1=a1+toc;
    
fprintf('Time needed to generate the linear systems = %0.1fs, for the solution = %.1fs\r',a1,a2)
fprintf([str,'\r'],err)



% Keep only the common coarse points (coarse structure)
for i=1:length(U)
    Pa2c{i}=Projector(grid{i}.GP,Coarse,LH);
%     P_bc{i}=ProjectorBC(Fine,grid{i}.GP,LH);
    Pf2a{i}=Projector(Fine,grid{i}.GP,LH);
%     Pf2a{i}=Projector_fine2anisotropic(grid{i}.GP,Fine,LH);
    Uc{i}=Pa2c{i}*U{i}(:);
end


ORD=[space_ord ord_time];


Coarse_12=[Coarse(1:2),Fine(3)];
Coarse_13=[Coarse(1),Fine(2),Coarse(3)];
Coarse_23=[Fine(1),Coarse(2:3)];
Pf2c{1}=Projector(Fine,Coarse_12,LH);
Pf2c{2}=Projector(Fine,Coarse_13,LH);
Pf2c{3}=Projector(Fine,Coarse_23,LH);




%intersection mesh 1 with mesh 2 without mesh 3
Pa2a{1,1}=Projector(grid{1}.GP,Coarse_12,LH);
Pa2a{1,2}=Projector(grid{2}.GP,Coarse_12,LH);
U_12{1}=Pa2a{1,1}*U{1}(:);
U_12{2}=Pa2a{1,2}*U{2}(:);
Nx=Coarse_12(1); Ny=Coarse_12(2); Nz=Coarse_12(3);
if TIME==1; Nt=Coarse_12(3); end
[var12,LH12,RH12,u12,~]=LS(Nx,Ny,Nz,Nt,TIME,DIM,stencil_Laplacian,stencil_First_ord);
[~,iAx12,b12,U12,grid12,~] = LinearSystemGenerator(var12,LH12,RH12,u12);
U_12{3}=iAx12(b12);
if TIME==1
    hh={h{1}(1:2),h{2}(1:2),grid12.h(1:2)};
else
    hh={h{1}(1:2),h{2}(1:2),grid12.h(1:2)};
end
Uextrap1=RichardsonExtrap(U_12,hh,LH,ORD);
Uextrap1=reshape(Uextrap1,Coarse_12);
% ERROR AFTER RICHARDSON EXTRAPOLATION
% err_x=norm(U_12{1}(:)-U12(:))/norm(U12(:))
% err_y=norm(U_12{2}(:)-U12(:))/norm(U12(:))
% err_extrap=norm(U12(:)-Uextrap1(:))/norm(U12(:))





%intersezione mesh 1 con mesh 3 senza i mesh 2
Pa2a{2,1}=Projector(grid{1}.GP,Coarse_13,LH);
Pa2a{2,2}=Projector(grid{3}.GP,Coarse_13,LH);
U_13{1}=Pa2a{2,1}*U{1}(:);
U_13{2}=Pa2a{2,2}*U{3}(:);
Nx=Coarse_13(1); Ny=Coarse_13(2); Nz=Coarse_13(3);
if TIME==1; Nt=Coarse_13(3); end
[var13,LH13,RH13,u13,~]=LS(Nx,Ny,Nz,Nt,TIME,DIM,stencil_Laplacian,stencil_First_ord);
[~,iAx13,b13,U13,grid13,~] = LinearSystemGenerator(var13,LH13,RH13,u13);
U_13{3}=iAx13(b13);
if TIME==1
    hh={[h{1}(1),dt{1}],[h{2}(1),dt{2}],[grid13.h(1),grid13.dt]};
else
    hh={h{1}([1 3]),h{3}([1 3]),grid13.h([1 3])};
end
Uextrap2=RichardsonExtrap(U_13,hh,LH,ORD);
Uextrap2=reshape(Uextrap2,Coarse_13);
% ERROR AFTER RICHARDSON EXTRAPOLATION
% err_x=norm(U_13{1}(:)-U13(:))/norm(U13(:))
% err_z=norm(U_13{2}(:)-U13(:))/norm(U13(:))
% err_extrap=norm(U13(:)-Uextrap2(:))/norm(U13(:))


ASD=Pf2c{2}*(Pf2c{2}'*Uextrap2(:)-Pf2c{1}'*(Pf2c{1}*(Pf2c{2}'*Uextrap2(:))));
% asd=Pf2c{2}'*Uextrap2(:);
asd=Pf2c{2}'*ASD;
asd2=Pf2c{1}*asd;
asd3=Pf2c{3}*asd;
if sum(asd2~=0)~=0 || sum(asd3~=0)~=0
    error('Error')
end




%intersezione mesh 2 con mesh 3 senza i mesh 1
Pa2a{3,1}=Projector(grid{2}.GP,Coarse_23,LH);
Pa2a{3,2}=Projector(grid{3}.GP,Coarse_23,LH);
U_23{1}=Pa2a{3,1}*U{2}(:);
U_23{2}=Pa2a{3,2}*U{3}(:);
Nx=Coarse_23(1); Ny=Coarse_23(2); Nz=Coarse_23(3);
if TIME==1; Nt=Coarse_23(3); end
[var23,LH23,RH23,u23,~]=LS(Nx,Ny,Nz,Nt,TIME,DIM,stencil_Laplacian,stencil_First_ord);
[~,iAx23,b23,U23,grid23,~] = LinearSystemGenerator(var23,LH23,RH23,u23);
U_23{3}=iAx23(b23);
if TIME==1
    hh={[h{2}(2),dt{2}],[h{3}(2),dt{3}],[grid23.h(2),grid23.dt]};
else
    hh={h{2}([2 3]),h{3}([2 3]),grid23.h([2 3])};
end
Uextrap3=RichardsonExtrap(U_23,hh,LH,ORD);
Uextrap3=reshape(Uextrap3,Coarse_23);
% ERROR AFTER RICHARDSON EXTRAPOLATION
% err_y=norm(U_23{1}(:)-U23(:))/norm(U23(:))
% err_z=norm(U_23{2}(:)-U23(:))/norm(U23(:))
% err_extrap=norm(U23(:)-Uextrap3(:))/norm(U23(:))


grid_sc={grid12,grid13,grid23};

CP1=Pa2c{1}*(Pf2a{1}*(Pf2c{1}'*Uextrap1(:)));
CP2=Pa2c{2}*(Pf2a{2}*(Pf2c{2}'*Uextrap2(:)));
CP3=Pa2c{3}*(Pf2a{3}*(Pf2c{3}'*Uextrap3(:)));
errr{1}=reshape(CP1-CP2,Coarse);
errr{2}=reshape(CP1-CP3,Coarse);

ord=5; ORDD={ord,ord,ord};
xx=[grid{end}.xx,grid{end}.tt];
for k=1:2
    e=padarray(errr{k},ones(1,length(grid{end}.GP)),0);
    if isempty(grid{end}.tt)==0
    	e=e(:,:,1:end-1);
    end
    p=spapi(ORDD,xx,e); %spline interp
    if k==1
        Uextrap2 = reshape(Uextrap2,grid_sc{k+1}.GP)+fnval(p,[grid_sc{k+1}.xx_intern,grid_sc{k+1}.tt_intern]);
    else
        Uextrap3 = reshape(Uextrap3,grid_sc{k+1}.GP)+fnval(p,[grid_sc{k+1}.xx_intern,grid_sc{k+1}.tt_intern]);
    end
end


Uextrap32=Pf2c{3}'*Uextrap3(:);
Proj1=Pf2c{1}'*(Pf2c{1}*Uextrap32);
ASD=Pf2c{3}*(Pf2c{3}'*Uextrap3(:)-Proj1);

Uextra=Pf2c{1}'*Uextrap1(:);
Uextra=Uextra+Pf2c{2}'*Uextrap2(:)-Pf2c{1}'*(Pf2c{1}*(Pf2c{2}'*Uextrap2(:)));
Uextra=Uextra+Pf2c{3}'*ASD;%Uextrap3(:);


% CheckIrreg1(Fine,Coarse,Uextra) % shows the strongest irregularitiy of the surface


%global check on the extrapolated solution:
Uex_fine=uex_fine(:); Uex_fine(Uextra==0)=0;
fprintf('Error after extrapolation : %0.3e\r',norm(Uextra-Uex_fine)/norm(Uex_fine))



% SKELETON 
[Ucorr,U2]=Correction3D(U,Pa2c,Pf2a,Pf2c,reshape(Uextra,Fine),grid,grid_sc); 
% CheckIrreg2(Fine,Coarse,Ucorr) % shows the strongest irregularity of the surface



POSS=Uextra==0;
Uex_fine=uex_fine(:); Uex_fine(Ucorr(:)==0)=0;
fprintf('Error SKELETON after extrapolation and correction = %0.2e\r',norm(Ucorr(:)-Uex_fine)/norm(Uex_fine))


% FILLING HOLES (UPSAMPLING)
tic
Ucorr=FillingHoles(Ucorr,Fine,Coarse,LH,A_fine,b_fine,TIME);

Ucorr=reshape(Ucorr,Fine);
e_up=norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:));
r_up=norm(b_fine(:)-A_fine*Ucorr(:))/norm(b_fine(:));
t_up=toc;
fprintf(['Error after UPsampling           = %0.1e,  Residual = %0.1e     (solution time = %.1fs)\r'],e_up,r_up,t_up) %norm(Uup(:)-Uf(:))/norm(Uf(:))


% ITERATE THE ALGORITHM
clear e
fig=0;



% corr=2; % correct: 1=x is more correct, 2=t is more correct,...
rich=0;
ORD=[space_ord ord_time];
for k=1:MAXIT
    res=b_fine(:)-A_fine*Ucorr(:);
    ress(k)=norm(res)/norm(b_fine(:));
    erro(k)=norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:));
    fprintf(['IT=%.0f - ERROR : %0.3e      RESIDUAL : %0.3e\r'],k-1,erro(end),ress(end)) %norm(Uup(:)-Uf(:))/norm(Uf(:))
    
    b_fine=b_fine(:); res=res(:);
    e_ex=0*res;%A_fine\res;
    if fig==1, figure(fig1); surf(X2,T2,reshape(res,Fine));  end
        
    e={};
    PP=Pf2a;
%     PP=P_bc;
    for proj=1:length(Pa2c)
        if proj==(length(Pa2c)-1) && TIME==1%proj==(length(P)-1) || proj<length(P)) && TIME==1
            r{proj}=PP{proj}*res;%/dt{1}*dt{2};   % P{proj}
        else
            r{proj}=PP{proj}*res;%*h{1}/h{2};   % P{proj}
        end
        
        e{proj}=A{proj}\r{proj};
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
        for proj=1:length(Pa2c)
            ec{proj}=Pa2c{proj}*e{proj};   % P{proj}
        end
        Uextrap=RichardsonExtrap(ec,h,dt,LH,ORD);
        Uextrap=reshape(Uextrap,Coarse);
        fprintf(['Accuracy of projected error (before extrapolation) : \r',str,'\r'],e2) %norm(Uup(:)-Uf(:))/norm(Uf(:))
        
        err_extrap=norm(Pf2a{end}*e_ex(:)-Uextrap(:))/norm(Pf2a{end}*e_ex(:));
        fprintf('Accuracy of COARSE error after extrapolation = %0.2e\r',err_extrap)
    else
        for ii=1:3
            E{ii}=Pf2a{ii}'*e{ii};
            E{ii}(POSS)=0;
        end
        Uextrap=E{1};
        
        %add E{2}
        for ii=2:3
            Err=reshape(Pa2c{1}*(Pf2a{1}*Uextrap),grid{end}.GP)-reshape(Pa2c{ii}*e{ii},grid{end}.GP);
            xx=[grid{end}.xx,grid{end}.tt]; %coarse mesh with boundaries
            ee=padarray(Err,ones(1,length(grid{end}.GP)),0);
            if isempty(grid{end}.tt)==0
                if length(size(U{1}))==2
                    ee=ee(:,1:end-1);
                elseif length(size(U{1}))==3
                    ee=ee(:,:,1:end-1);
                elseif length(size(U{1}))==4
                    ee=ee(:,:,:,1:end-1);
                elseif length(size(U{1}))==5
                    ee=ee(:,:,:,:,1:end-1);
                end
            end
            p=spapi(ORDD,xx,ee); %spline interp
            e{ii}=e{ii}+reshape(fnval(p,[grid{ii}.xx_intern,grid{ii}.tt_intern]),[],1);  
%             CheckIrreg1(grid{ii}.GP,Coarse,e{ii})
        end
        p1=abs(Uextrap)<10^-15;       c1=Pf2a{2}'*e{2};
        Uextrap(p1 & POSS==0)=c1(p1 & POSS==0);
    end
    
    % Correction 
    [Ecorr,~]=Correction3D(e,Pa2c,Pf2a,Pf2c,reshape(Uextrap,Fine),grid,grid_sc); 
%     CheckIrreg2(Fine,Coarse,Ecorr)  %shows the strongest irregularity of the surface


    % Filling holes
    Ecorr=FillingHoles(Ecorr,Fine,Coarse,LH,A_fine,zeros(size(b_fine)),TIME);
    

    %Apply correction
    f=@(c)norm(b_fine(:)-A_fine*reshape(Ucorr(:)+c*Ecorr(:),[],1));
    OPTIONS.TolX=10^-10;
    c2=fminbnd(f,0,10,OPTIONS);
    Ucorr=Ucorr(:)+c2*Ecorr(:);
    e{end+1}=c2*Ecorr;    
end
fprintf(['IT=%.0f - ERROR : %0.3e      RESIDUAL : %0.3e\r'],k,norm(Ucorr(:)-uex_fine(:))/norm(uex_fine(:)),norm(res)/norm(b_fine(:))) %norm(Uup(:)-Uf(:))/norm(Uf(:))




return














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

if 1==2
    [blocks,FLAG,RELRES,ITER,RESVEC]=gmres(A2,f2-P_f2holes*A_fine*Ucorr(:),100,10^-10,[],@(r)reshape(B1\reshape(r,hole(1),[]),[],1));
else
    % This is slow and can be largely improved (enough for testing)
    b=f2-P_f2holes*A_fine*Ucorr(:);
    blocks=zeros(size(b));
    
    if TIME==1
        B1=A2(1:prod(hole),1:prod(hole));
        for it=1:2
            blocks=blocks+reshape(B1\reshape(b-A2*blocks,prod(hole),[]),[],1);
        end
        [blocks,FLAG,RELRES,ITER,RESVEC]=gmres(A2,f2-P_f2holes*A_fine*Ucorr(:),100,10^-10,[],@(r)reshape(B1\reshape(r,prod(hole),[]),[],1),[],blocks);
    else
        v=hole(1)*prod(space(2:end));
        v=prod(space(2:end));
        B1=A2(1:v,1:v);
        w=1;
        for it=1:20
            blocks=blocks+w*reshape(B1\reshape(b-A2*blocks,v,[]),[],1);
            if it>10
                w=rand;
            end
        end
        [blocks,FLAG,RELRES,ITER,RESVEC]=gmres(A2,b,100,10^-10,[],@(r)reshape(B1\reshape(r,v,[]),[],1),[],blocks);
    end
end
Ucorr=Ucorr(:)+P_f2holes'*blocks;
end





function P=Projector_Fine2Holes(Fine,Coarse,LH)
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




function [Ucorr,U]=Correction3D(U,Pa2cc,Pf2a,Pf2c,Uextrap,grid,grid_sc)
% First adjust the tree surfaces s.t. the very coarse points are corrected
for i=1:length(U)
    Err{i}=reshape(Pa2cc{i}*(Pf2a{i}*Uextrap(:)),grid{end}.GP)-reshape(Pa2cc{i}*U{i}(:),grid{end}.GP);
end

ord=5; %spline interpolation order
for kk=1:length(grid{end}.xx)
    ORD{kk}=ord;
end
if isempty(grid{end}.tt)==0
    ORD{kk+1}=ord;
end

for i=1:length(U)-1 
    for j=1:length(Pf2c)
        if (i==1 && (j==1 || j==2)) || (i==2 && (j==1 || j==3)) || (i==3 && (j==2 || j==3))
            temp=Pf2c{j}*(Pf2a{i}'*U{i}(:));
            temp1=Pf2c{j}*Uextrap(:);
            Err{j}=reshape(temp1-temp,grid_sc{j}.GP);
            xx=[grid_sc{j}.xx,grid_sc{j}.tt]; %coarse mesh with boundaries
            e=padarray(Err{j},ones(1,length(grid_sc{j}.GP)),0);
            if isempty(grid{end}.tt)==0
                if length(size(e))==2
                    e=e(:,1:end-1);
                elseif length(size(e))==3
                    e=e(:,:,1:end-1);
                elseif length(size(e))==4
                    e=e(:,:,:,1:end-1);
                elseif length(size(e))==5
                    e=e(:,:,:,:,1:end-1);
                end
            end
            p=spapi(ORD,xx,e); %spline interp
            U{i} = U{i}(:)+reshape(fnval(p,[grid{i}.xx_intern,grid{i}.tt_intern]),[],1);
        end
    end
end

for i=1:length(U)-1
    if i>1
        pos=Ucorr==0;
        update=Pf2a{i}'*U{i}(:);
        Ucorr(pos)=update(pos);
    else
        Ucorr=Pf2a{i}'*U{i}(:);
    end
end
end


% profile viewer
function Uextrap=RichardsonExtrap(U,h,LH,ORD)
% Removes lowest order terms
for i=1:length(U)-1
    for j=1:length(h{i})
        LS(j,i)=h{i}(j)^ORD(1);% generating linear system
    end
    if i<=length(U)-1
        b(i,1)=h{end}(i)^ORD(1);
    end
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
%                 PP=speye(Fine(i));
                s=Fine(i)/Coarse(i);
%                 if mod(s,2)==0 % s is integer
%                     PP=PP(s/2:s:Fine(end),:);
%                 else
%                     PP=PP((s+1)/2:s:Fine(end),:);
%                 end
                PP=speye(Fine(i));
                for j=1:Coarse(i)
                    if j==1
                        PP=sparse(1,Fine(i));
                        PP(j,s)=1;
                        PP(j,1)=1;
%                     elseif j==Coarse(i)
%                         PP(j,s*j)=1;
%                         PP(j,s*j)=1;
                    else
                        PP(j,s*j)=1;
                        PP(j,s*(j-1))=1;
                    end
                end
                PP=PP/2;
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




function [var,LH,RH,u,ord_time]=LS(Nx,Ny,Nz,Nt,TIME,DIM,stencil_Laplacian,stencil_First_ord)
    global discr nu1 nu2 nu3 d1 d2 d3 uu 

    if TIME==1
        if DIM==1
            var={'x',Nx;'t',Nt};
            u=char(uu);
            f=char(diff(uu,1,'t')-(d1.*diff(uu,2,'x')-nu1.*diff(uu,1,'x')));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                [],f,f};
        elseif DIM==2
            var={'x',Nx;'y',Ny;'t',Nt};
            u=char(uu);
            f=char(diff(uu,1,'t')-(d1.*diff(uu,2,'x')-nu1.*diff(uu,1,'x')+d2.*diff(uu,2,'y')-nu2.*diff(uu,1,'y')));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'y',char(d2),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                'y',char(nu2),@(h)stencil_First_ord/h;
                [],f,f};
        end
    else
        if DIM==2
            var={'x',Nx;'y',Ny};
            u=char(uu);
            f=char(-d1.*diff(uu,2,'x')-d2.*diff(uu,2,'y')+nu1.*diff(uu,1,'x')+nu2.*diff(uu,1,'y'));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'y',char(d2),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                'y',char(nu2),@(h)stencil_First_ord/h;
                [],f,f};    
        elseif DIM==3
            var={'x',Nx;'y',Ny;'z',Nz};
            u=char(uu);
            f=char(-d1.*diff(uu,2,'x')-d2.*diff(uu,2,'y')-d3.*diff(uu,2,'z')+nu1.*diff(uu,1,'x')+nu2.*diff(uu,1,'y')+nu3.*diff(uu,1,'z'));
            RH={'x',char(d1),@(h)stencil_Laplacian/h^2;
                'y',char(d2),@(h)stencil_Laplacian/h^2;
                'z',char(d3),@(h)stencil_Laplacian/h^2;
                'x',char(nu1),@(h)stencil_First_ord/h;
                'y',char(nu2),@(h)stencil_First_ord/h;
                'z',char(nu3),@(h)stencil_First_ord/h;
                [],f,f};
        end
    end
    if TIME==1
        if discr==1
            LH={'CN'}; ord_time=2;
        elseif discr==2
            LH={'IE'}; ord_time=1;
        else
            LH={'EE'}; ord_time=1;
        end
    else
        ord_time=0;
        LH={};
    end


end



function CheckIrreg1(Fine,Coarse,U)
U=reshape(U,Fine);
n=ones(size(Fine));
for i=1:length(Fine)
    while n(i)*(Coarse(i)+1)-1~=Fine(i) && n(i)<10^5
        n(i)=n(i)+1;
    end
end
if any(n==10^5-1)
    error('something went wrong')
end

v1=[];
for i=1:Coarse(1)-1
for j=1:Coarse(2)-1
    d(i,j)=max(diff(reshape(U(i*n(1),j*n(2),:),[],1)));
    v1(end+1,1)=i;
    v1(end,2)=j;
end
end
d=abs(d);
[pos1,val(1)]=find(d==max(d(:)));
v(1,1:2)=v1(min(pos1),:);

v2=[];
for i=1:Coarse(1)-1
for j=1:Coarse(3)-1
    d(i,j)=max(diff(reshape(U(i*n(1),:,j*n(3)),[],1)));
    v2(end+1,1)=i;
    v2(end,2)=j;
end
end
d=abs(d);
[pos2,val(2)]=find(d==max(d(:)));
v(2,1:2)=v2(min(pos2),:);

v3=[];
for i=1:Coarse(2)-1
for j=1:Coarse(3)-1
    d(i,j)=max(diff(reshape(U(:,i*n(2),j*n(3)),[],1)));
    v3(end+1,1)=i;
    v3(end,2)=j;
end
end
d=abs(d);
[pos3,val(3)]=find(d==max(d(:)));
v(3,1:2)=v3(min(pos3),:);

pos=find(val==max(val));
pos=min(pos);
v=v(pos,:);
i=v(1);
j=v(2);

figure
if pos==1
    plot(reshape(U(i*n(1),j*n(2),:),[],1))
elseif pos==2
    plot(reshape(U(i*n(1),:,j*n(3)),[],1))
else
    plot(reshape(U(:,i*n(2),j*n(3)),[],1))
end

end





function CheckIrreg2(Fine,Coarse,U)
U=reshape(U,Fine);
n=ones(size(Fine));
for i=1:length(Fine)
    while n(i)*(Coarse(i)+1)-1~=Fine(i) && n(i)<10^5
        n(i)=n(i)+1;
    end
end
if any(n==10^5-1)
    error('something went wrong')
end


v1=[]; c=0;
for k=1:6
if k==1
    c1=Coarse(1)-1;
    f1=Fine(2);
elseif k==2
    c1=Coarse(2)-1;
    f1=Fine(1);
elseif k==3 
    c1=Coarse(1)-1;
    f1=Fine(3);
elseif k==4
    c1=Coarse(3)-1;
    f1=Fine(1);
elseif k==5
    c1=Coarse(2)-1;
    f1=Fine(3);
else
    c1=Coarse(3)-1;
    f1=Fine(2);
end
for i=1:c1
for j=1:f1
    c=c+1;
    if k==1
        d(c)=max(diff(reshape(U(i*n(1),j,:),[],1)));
        v1(c,1)=i*n(1);
        v1(c,2)=j;
        v1(c,3)=0;
    elseif k==2
        d(c)=max(diff(reshape(U(j,i*n(2),:),[],1)));
        v1(c,1)=j;
        v1(c,2)=i*n(2);
        v1(c,3)=0;
    elseif k==3
        d(c)=max(diff(reshape(U(i*n(1),:,j),[],1)));
        v1(c,1)=i*n(1);
        v1(c,2)=0;
        v1(c,3)=j;
    elseif k==4
        d(c)=max(diff(reshape(U(j,:,i*n(3)),[],1)));
        v1(c,1)=j;
        v1(c,2)=0;
        v1(c,3)=i*n(3);
    elseif k==5
        d(c)=max(diff(reshape(U(:,i*n(2),j),[],1)));
        v1(c,1)=0;
        v1(c,2)=i*n(2);
        v1(c,3)=j;
    elseif k==6
        d(c)=max(diff(reshape(U(:,j,i*n(3)),[],1)));
        v1(c,1)=0;
        v1(c,2)=j;
        v1(c,3)=i*n(3);
    end
end
end
end

pos=find(d==max(d),1);
% pos=min(pos);
v1=v1(pos,:);
figure
if v1(1)==0
    plot(reshape(U(:,v1(2),v1(3)),[],1))
elseif v1(2)==0
    plot(reshape(U(v1(1),:,v1(3)),[],1))
else
    plot(reshape(U(v1(1),v1(2),:),[],1))
end
end

% figure
% for i=1:79
%     surf(reshape(U(i,:,:),79,79))
%     pause(.3)
% end



