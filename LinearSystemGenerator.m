function [Ax,iAx,b,uex,grid,A] = LinearSystemGenerator(vars,LH,RH,u)
% Generator of the linear system from the discretization of the linear PDE
% du/dt=F(u,x,y,t,...); or 0=F(u,x,y,...);
%
% DISCRETIZATION OVER UNIFORM MESHES WITH N=[N_1,N_2,...] GRID POINTS OVER
% THE MESHES IN EACH DIRECTION var{i,1}.
% var = contains the label of the coordinate axis and the amount of grid
%       poitns in that direction, e.g.:
%       var={'x',10;'y',20;'t',100} -> 10 points on the mesh in x, 20 in y and 100 in time
% LH = left  hand side is devoted to time du/dt 
% e.g. - LH={'CN'}: time-dependent PDE, discretized through Crank Nicolson
%      - LH={}: stationary PDE, no time dependency
% RH = right hand side is devoted to the rest of the pde, e.g. d^2u/dx^2 +...
% e.g. - RH{'x',c1,sten1 ;  -> some derivative with respect to "x" where
%                              - "c" is the constant/function multiplying the
%                                derivative (char type e.g. '1+x.*t')
%                              - "sten" is the stencil, e.g. sten=@(h)[-1 2 -1]/h^2
%                                for the laplacian -d^2 u/dx^2
%           'x',c2,sten2;   -> some other derivative with respect to x
%           'y',c3,sten3}   -> some other derivative with respect to y
% 
% WARNING : Mixed derivatives not implemented yet
%
% ASSUMPTIONS:
% - Domain is considered to be the square [0 1]^d
%

%CHECK
if size(RH,2)~=3
    error("RH must be a cellarray with 3 columns, e.g. RH={'x',c1,sten1;'y',c3,sten3}")
end
if length(LH)~=1 && ~isempty(LH)
    error("LH must be a cellarray with 1 element, e.g. RH={'CN'} or IE,...")
elseif isempty(LH)==0 && strcmp(LH,'IE')==0 && strcmp(LH,'EE')==0 && strcmp(LH,'CN')==0
    error('Only IE (implicit euler), EE (explicit euler) and CN (crank nicolson) are implemented')
end
if size(vars,2)~=2
error("'var' must have 2 columns, the first contains the label of the coordinate axes and the second the amount of grid point the respective direction")
end

vars=sortrows(vars);
time=0; space_vars=vars;
if any(strcmp(vars(:),'t'))
    time=1;
    time_var=space_vars(strcmp(space_vars(:,1),'t'),:);
    space_vars(strcmp(space_vars(:,1),'t'),:)=[];
end
if time
    dim=length(vars)-1;
else
    dim=length(vars);
end

%looking for forcing term (independent of u)
counter=0;
for i=1:size(RH,1)
    if strcmp(RH{i,1},'') || isempty(RH{i,1}) % this part is not a derivative
        counter=counter+1;
        f=AddDots(RH{i,2});
        if ischar(f)==0
            error('Every function has to be in char type')
        end
    end
end
if counter>1
    error('There should be only one forcing term')
end


%generate steps
for i=1:dim
    GP(i)=space_vars{i,2}; %grid points
    GP2(i)=GP(i)+2; 
    h(i)=1/(GP(i)+1);
end
if time
    GP(dim+1)=time_var{1,2};
    GP2(dim+1)=GP(dim+1)+1;
    dt=1/GP(dim+1);
else
    dt=[];
end
grid.GP=GP;
grid.h=h;
grid.dt=dt;

%generate meshes
str=[]; str2=[];
for i=1:dim
    xx{i}=0:h(i):1;
    xx_intern{i}=h(i):h(i):(1-h(i));
    if i==1
        str=[str,'xx{',num2str(i),'}'];
        str2=[str2,'X{',num2str(i),'}'];
    else
        str=[str,',xx{',num2str(i),'}'];
        str2=[str2,',X{',num2str(i),'}'];
    end
end
grid.xx=xx;
grid.xx_intern=xx_intern;

tt=[];tt_intern=[]; T=[];
if time
    tt=0:dt:1;
    tt_intern=dt:dt:1;
    eval(['[',str2,',T]= ndgrid(',str,',tt);'])
else
    eval(['[',str2,']= ndgrid(',str,');'])
end
grid.tt=tt;
grid.tt_intern=tt_intern;
grid.X=X;
grid.T=T;

% solution
u=AddDots(u);
uex=Feval(u,X,T,space_vars,0); %variable coefficient multiplying the differential operator


% Building the matrices in space
for i=1:dim 
    k=0; j=0;
    for jj=1:size(RH,1)
        if strcmp(space_vars{i,1},RH{jj,1})
            j=j+1;
            k=k+1;
            if time==1
                t=T;
            end
            c=Feval(AddDots(RH{jj,2}),X,T,space_vars,1); %variable coefficient multiplying the differential operator
            c2=Feval(AddDots(RH{jj,2}),X,T,space_vars,0); %variable coefficient multiplying the differential operator
            c=c(:); c2=c2(:);
            if all(abs(diff(c))<10^-14)
                %             |-------------      Toeplitz operator   ----------------------|
                Adim{i,j,1}=spdiags(ones(space_vars{i,2},1)*RH{jj,3}(h(i)),-1:1,space_vars{i,2},space_vars{i,2}); %storing the operator
                Adim{i,j,2}=sparse(c(1)); %storing the coefficient which multiplies the operator
                Adim2{i,j,1}=spdiags(ones(space_vars{i,2}+2,1)*RH{jj,3}(h(i)),(-1:1),space_vars{i,2}+2,space_vars{i,2}+2); %storing the larger operator.. needed for BC
                Adim2{i,j,2}=sparse(c2(1)); %storing the coefficient which multiplies the operator
            else
                %           |-------------      Toeplitz operator   ----------------------|
                Adim{i,j,1}=spdiags(ones(space_vars{i,2},1)*RH{jj,3}(h(i)),-1:1,space_vars{i,2},space_vars{i,2}); %storing the operator
                Adim{i,j,2}=spdiags(c(:),0,prod([vars{:,2}]),prod([vars{:,2}])); %storing the coefficient which multiplies the operator
                Adim2{i,j,1}=spdiags(ones(space_vars{i,2}+2,1)*RH{jj,3}(h(i)),(-1:1),space_vars{i,2}+2,space_vars{i,2}+2); %storing the larger operator.. needed for BC
                Adim2{i,j,2}=spdiags(c2(:),0,prod(GP2),prod(GP2)); %storing the coefficient which multiplies the operator
            end
        end
    end
    if time
        if strcmp(LH,'IE') || strcmp(LH,'EE')
            At=spdiags(ones(time_var{1,2},1)*[-1 1]/dt,-1:0,time_var{1,2},time_var{1,2});
            At2=spdiags(ones(time_var{1,2}+1,1)*[-1 1]/dt,-1:0,time_var{1,2}+1,time_var{1,2}+1);
        elseif strcmp(LH,'CN')
            At=spdiags(ones(time_var{1,2},1)*[-1 1]/(dt),-1:0,time_var{1,2},time_var{1,2});
            At2=spdiags(ones(time_var{1,2}+1,1)*[-1 1]/(dt),-1:0,time_var{1,2}+1,time_var{1,2}+1);
        else
            error('Not implemented')
        end
    else
        At=[];
        At2=[];
    end
end


% Boundary conditions
for i=1:dim+time
    BC{i}=[];
    for j=1:dim
        BC{i}=[BC{i},':,'];
    end
    if time
        BC{i}=[BC{i},':,'];
    end
    BC{i}=BC{i}(1:end-1); %removing final comma
end
BC0=BC; BC1=BC; clear BC
for i=1:dim+time
    BC0{i}(2*i-1)='1';
    BC1{i}=[BC1{i}(1:2*i-2),'end',BC1{i}(2*i:end)];
end

if isempty(LH)
    b=Feval(f,X,T,space_vars,1);
elseif strcmp(LH,'IE')
    b=Feval(f,X,T,space_vars,1);
elseif strcmp(LH,'EE')
    b=Feval(f,X,T-dt,space_vars,1);
elseif strcmp(LH,'CN')
    b=Feval(f,X,T-dt/2,space_vars,1);
end
grid.f=b;

boundaries=[];
for i=1:size(space_vars,1)
    boundaries=[boundaries,'2:end-1,'];
end
if time
    boundaries=[boundaries,'2:end,'];
end
eval(['ubc=uex;ubc(',boundaries(1:end-1),')=0;'])


fubc=A_times_x(Adim2,At2,ubc,GP2,LH);
fubc=reshape(fubc,GP2);
eval(['fubc=fubc(',boundaries(1:end-1),');'])
eval(['fubc(',boundaries(1:end-1),')=0;'])

b=b-fubc;





% eval(['b(',BC0{end},')=b(',BC0{end},')+b0(:);'])


% if time % condition at t0=0
%     str=[];
%     for i=1:length(space_vars)
%         str=[str,space_vars{i},'=X{',num2str(i),'}(',BC0{end},')'];
%     end
%     
%     boundaries=[];
%     for i=1:size(space_vars,1)
%         boundaries=[boundaries,'2:end-1,'];
%     end
%     boundaries=[boundaries,'1'];
%     eval(['u0=uex(',boundaries,');'])
%     if strcmp(LH,'IE')
%         t=0;
%         b0=Feval(f,X,0,space_vars,1);
%     elseif strcmp(LH,'EE')
%         
%     elseif strcmp(LH,'CN')
%         
%     else
%         error('Not implemented')
%     end
% %     for i=1:size(space_vars,1)
% %         boundaries=[boundaries,'2:end-1,'];
% %     end
% %     eval(['ubc=uex;ubc(',,')=0;'])
% %     eval(['b(',BC0{end},')=b(',BC0{end},')+b0(:);'])
% end
% for i=1:dim
%     eval(['b(',BC0{i},')=;'])
%     eval(['b(',BC1{i},')=;'])
% end










uex=Feval(u,X,T,space_vars,1); %variable coefficient multiplying the differential operator
Ax=@(x)A_times_x(Adim,At,x,GP,LH);
if nargout==6
%     tic
    A=Assemble(Adim,At,GP,LH);
%     toc
%     'asd'
else
    A=[];
end
iAx=@(x)inverseA_times_x(Adim,At,x,GP,LH,A);



% TEST
check=0;
if check
    fx=Feval(f,X,T,space_vars,1);
    fx2=Feval(f,X,T-dt/2,space_vars,1);
    tic
    f2=A_times_x(Adim,At,uex,GP,LH);
    toc
    tic
    f3=Assemble(Adim,At,GP,LH)*uex(:);
    toc
    norm(fx(:)-f2(:))/norm(fx(:))
    norm(fx(:)-f3(:))/norm(f2(:))
    norm(fx2(:)-f3(:))/norm(f2(:))
end

end




function fx=Feval(f,X,T,space_vars,inside)
    if ischar(f)
        str=[];
        if inside % removing boundaries
            boundaries=[];
            for i=1:length(X)
                boundaries=[boundaries,'2:end-1,'];
            end
            if isempty(T)==0
                boundaries=[boundaries,'2:end,'];
            end
            for i=1:length(X)
                str=[str,space_vars{i},'=X{',num2str(i),'}(',boundaries(1:end-1),');'];
            end
            if isempty(T)==0
                str=[str,'t=T(',boundaries(1:end-1),');'];
            end
        else
            for i=1:length(X)
                str=[str,space_vars{i},'=X{',num2str(i),'};'];
            end
            if isempty(T)==0
                str=[str,'t=T;'];
            end
        end
        eval([str,'fx=',f,';'])
    else
        error('"f" has to be a char')
    end
end




function f=AddDots(f)
if ischar(f)==0
    error('f must be char')
end
for i=length(f):-1:1
    if strcmp(f(i),'*') || strcmp(f(i),'/') || strcmp(f(i),'^')
        if strcmp(f(i-1),'.')==0 
            f=[f(1:i-1),'.',f(i:end)];
        end
    end
end
end




function y=A_times_x(Adim,At,x,GP,LH)
y=zeros(prod(GP),1);
x=reshape(x,GP);
if isempty(LH)
    shift=1; %blocks on the main diagonal of the all-at-once linear system
elseif strcmp(LH,'IE')
    shift=1; %blocks on the main diagonal of the all-at-once linear system
elseif strcmp(LH,'EE')
    shift=speye(GP(1));
    for ii=2:length(GP)-1
        shift=kron(speye(GP(ii)),shift);
    end
    shift=kron(spdiags(ones(GP(end),1),-1,GP(end),GP(end)),shift);%blocks on the sub-diagonal of the all-at-once linear system
elseif strcmp(LH,'CN')
    shift=speye(GP(1));
    for ii=2:length(GP)-1
        shift=kron(speye(GP(ii)),shift);
    end
    shift=kron(spdiags(ones(GP(end),2)/2,-1:0,GP(end),GP(end)),shift);%blocks on both the sub and main diagonals of the all-at-once linear system
end   
for i=1:size(Adim,1)
    %*spdiags(ones(GP(end),1),-1,GP(end),GP(end))'
    for j=1:size(Adim,2)
        if ~isempty(Adim{i,j,1})
            if length(GP)==2
                if i==1
                    y=y-Adim{i,j,2}*reshape(Adim{i,j,1}*x,[],1); %blocks on the main diagonal of the all-at-once linear system
                else
                    y=y-Adim{i,j,2}*reshape(x*Adim{i,j,1}',[],1); %blocks on the main diagonal of the all-at-once linear system
                end
            elseif length(GP)==3
                y2=zeros(GP);
                if i==1
                    for ii=1:GP(3)
                        y2(1:GP(1),1:GP(2),ii)=Adim{i,j,1}*reshape(x(:,:,ii),[GP(1:2)]);
                    end
                elseif i==2
                    for ii=1:GP(1)
                        y2(ii,1:GP(2),1:GP(3))=Adim{i,j,1}*reshape(x(ii,:,:),[GP(2:end)]);
                    end
                else
                    for ii=1:GP(1)
                        y2(ii,1:GP(2),1:GP(3))=reshape(x(ii,:,:),[GP(2:end)])*Adim{i,j,1}';
                    end
                end
                y=y-Adim{i,j,2}*y2(:);
            elseif length(GP)==4
                if i==1
                    for ii=1:GP(3), for jj=1:GP(4)
                        y2(1:GP(1),1:GP(2),ii,jj)=Adim{i,j,1}*reshape(x(:,:,ii,jj),[GP(1:2)]);
                    end,end
                elseif i==2
                    for ii=1:GP(3), for jj=1:GP(4)
                        y2(1:GP(1),1:GP(2),ii,jj)=reshape(x(:,:,ii,jj),[GP(1:2)])*Adim{i,j,1}';
                    end,end
                elseif i==3
                    for ii=1:GP(1), for jj=1:GP(4)
                        y2(ii,1:GP(2),1:GP(3),jj)=reshape(x(ii,:,:,jj),[GP(2:3)])*Adim{i,j,1}';
                    end,end
                else
                    for ii=1:GP(1), for jj=1:GP(2)
                        y2(ii,jj,1:GP(3),1:GP(4))=reshape(x(ii,jj,:,:),[GP(3:4)])*Adim{i,j,1}';
                    end,end
                end
                y=y-Adim{i,j,2}*y2(:);
            else
                %define the computation according to the properties of the tensor product
                error('Not implemented yet')
            end
        end
    end
end
% subplot(1,2,1)
% surf(reshape(y,GP(1),GP(2)))
y=shift*y;
% subplot(1,2,2)
% surf(reshape(y,GP(1),GP(2)))
if ~isempty(At)
    if length(GP)==2
        y2=reshape(x*At',[],1);
    elseif length(GP)==3
        y2=zeros(GP);
        for i=1:GP(1)
            y2(i,1:GP(2),1:GP(3))=reshape(x(i,:,:),[GP(2:end)])*At';
        end
    elseif length(GP)==4
        y2=zeros(GP);
        for i=1:GP(1)
        for j=1:GP(2)
            y2(i,j,1:GP(3),1:GP(4))=reshape(x(i,j,:,:),[GP(3:end)])*At';
        end
        end
    else
        %define the computation according to the properties of the tensor product
        error('Not implemented yet')
    end
    y=y+y2(:);
end

end


function iAx=inverseA_times_x(Adim,At,x,GP,LH,A)
x=x(:);
if exist('A','var')==0
    A=[];
end   
if isempty(At) || numel(At)==1 %no time
    if isempty(A)
        A=Assemble(Adim,At,GP,LH);
    end
    iAx=A\x; %preconditioned iterative methods can be used here
%     tic
%     [A,B]=Assemble(Adim,At,GP,LH);
%     iAx=gmres(A,x,100,10^-7,[],spdiags(diag(B),0,length(B),length(B)));
%     toc
%     'asd'
else
    if isempty(A)
        %time stepping exact solver
        n=prod(GP(1:end-1));
        for k=1:GP(end)
            if k==2
                counter=0;
                for i=1:size(Adim,1)
                for j=1:size(Adim,2)
                    if isempty(Adim{i,j,2})
                        counter=counter+1; 
                    elseif numel(Adim{i,j,2})==1
                        counter=counter+1; 
                    end
                end
                end
                cond=counter==(size(Adim,1)*size(Adim,2));
            end
            if k==1
    %             [A,Asub]=Blocks(Adim,At,GP,LH,k);
                A=Assemble(Adim,At,GP,LH);
            end
            index=((k-1)*n+1):(k*n);
            index2=index-n;
            if k==1
                iAx(index,1)=A(index,index)\x(index);
            else
                iAx(index,1)=A(index,index)\(x(index)-A(index,index2)*iAx(index2,1));
            end
        end
    else
        %time stepping exact solver
        n=prod(GP(1:end-1));
        for k=1:GP(end)
            index=((k-1)*n+1):(k*n);
            index2=index-n;
            if k==1
                iAx(index,1)=A(index,index)\x(index);
            else
                iAx(index,1)=A(index,index)\(x(index)-A(index,index2)*iAx(index2,1));
            end
        end
    end
end
%     if ischar(f)
%         str=[];
%         for i=1:length(X)
%             str=[str,space_var{i},'=X{',num2str(i),'};'];
%         end
%         if isempty(T)==0
%             str=[str,'t=T;'];
%         end
%         eval([';fx=',f,';'])
%         fx=fx(:);
%     else
%         error('')
%     end
end


function [A,Asub]=Blocks(Adim,At,GP,LH,t)
% Compute the diagonal and subdiagonal blocks at time t
if strcmp(LH,'IE')
    compute=[1 0];
elseif strcmp(LH,'EE')
    compute=[0 1];
elseif strcmp(LH,'CN')
    compute=[1 1];
end

n=prod(GP(1:end-1));
index=((At-1)*n+1):(At*n);

time=(isempty(At)==0 && numel(At)>1);
if time==0
    error('The matrix should depend on time')
end
time=0;
A=sparse(prod(GP(1:end-1)),prod(GP(1:end-1)));
Asub=sparse(prod(GP(1:end-1)),prod(GP(1:end-1)));
for i=1:size(Adim,1)
    for j=1:size(Adim,2)
        if isempty(Adim{i,j,1})==0
            if (size(Adim,1)+time)==1
                A2=Adim{i,j,1};
            elseif (size(Adim,1)+time)==2
                if i==1
                    A2=kron(speye(GP(2)),Adim{i,j,1});
%                     A2=kron(shift,Adim{i,j,1});
                else
                    A2=kron(Adim{i,j,1},speye(GP(1)));
                end
            elseif (size(Adim,1)+time)==3
                if i==1
                    A2=kron(speye(GP(3)),kron(speye(GP(2)),Adim{i,j,1}));
                elseif i==2
                    A2=kron(speye(GP(3)),kron(Adim{i,j,1},speye(GP(1))));
                else
                    A2=kron(Adim{i,j,1},kron(speye(GP(2)),speye(GP(1))));
                end
            elseif (size(Adim,1)+time)==4
                if i==1
                    A2=kron(kron(speye(GP(4)),kron(speye(GP(3)),speye(GP(2)))),Adim{i,j,1});
                elseif i==2
                    A2=kron(kron(speye(GP(4)),kron(speye(GP(3)),Adim{i,j,1})),speye(GP(1)));
                elseif i==3
                    A2=kron(kron(speye(GP(4)),kron(Adim{i,j,1},speye(GP(2)))),speye(GP(1)));
                else
                    A2=kron(kron(Adim{i,j,1},kron(speye(GP(3)),speye(GP(2)))),speye(GP(1)));
                end
            else
                error('Not implemented yet')
            end
            
            if t==1 
                if numel(Adim{i,j,2})>1
                    if compute(1)==1
                        A=A+Adim{i,j,2}(index,index)*A2;
                    end
                    if compute(2)==1
                        Asub=sparse(0);
                    end
                else
                    if compute(1)==1
                        A=A+Adim{i,j,2}*A2;
                    end
                    if compute(2)==1
                        Asub=sparse(0);
                    end
                end
            else
                if numel(Adim{i,j,2})>1
                    if compute(1)==1
                        A=A+Adim{i,j,2}(index,index)*A2;
                    end
                    if compute(2)==1
                        Asub=Asub+Adim{i,j,2}(index-n,index-n)*A2;
                    end
               	else
                    if compute(1)==1
                        A=A+Adim{i,j,2}*A2;
                    end
                    if compute(2)==1
                        Asub=Asub+Adim{i,j,2}*A2;
                    end
                end
            end
        end
%         'kron(kron(kron(a,b),c),d)'
    end    
end
if strcmp(LH,'IE')
    A=spdiags(diag(At),0,length(A),length(A))-A;
    Asub=At(2,1);%spdiags(diag(At,-1),0,length(At),length(At));
elseif strcmp(LH,'EE')
    A=At(1);%spdiags(diag(At),0,length(At),length(At));
    Asub=spdiags(diag(At,-1),-1,length(Asub),length(Asub))-Asub;
elseif strcmp(LH,'CN')
    A=spdiags(diag(At),0,length(A),length(A))-A/2;
    Asub=spdiags([diag(At,-1)],-1,length(Asub),length(Asub))-Asub/2;
end
if t==1
    Asub=[];
end

end




function [A,B]=Assemble(Adim,At,GP,LH)
if isempty(LH) || isempty(At)|| numel(At)==1
    shift1=1; %blocks on the main diagonal of the all-at-once linear system
elseif strcmp(LH,'IE')
    shift1=1; %blocks on the main diagonal of the all-at-once linear system
elseif strcmp(LH,'EE')
    shift1=speye(GP(1));
    for ii=2:length(GP)-1
        shift1=kron(speye(GP(ii)),shift1);
    end
    shift1=kron(spdiags(ones(GP(end),1),-1,GP(end),GP(end)),shift1);%blocks on the sub-diagonal of the all-at-once linear system
elseif strcmp(LH,'CN')
    shift1=speye(GP(1));
    for ii=2:length(GP)-1
        shift1=kron(speye(GP(ii)),shift1);
    end
    shift1=kron(spdiags(ones(GP(end),2)/2,-1:0,GP(end),GP(end)),shift1);%blocks on both the sub and main diagonals of the all-at-once linear system
end    
   
time=(isempty(At)==0 && numel(At)>1);
A=sparse(prod(GP),prod(GP));
for i=1:size(Adim,1)
    for j=1:size(Adim,2)
        if isempty(Adim{i,j,1})==0
            if (size(Adim,1)+time)==1
                A2=Adim{i,j,1};
            elseif (size(Adim,1)+time)==2
                if i==1
                    A2=kron(speye(GP(2)),Adim{i,j,1});
%                     A2=kron(shift,Adim{i,j,1});
                else
                    A2=kron(Adim{i,j,1},speye(GP(1)));
                end
            elseif (size(Adim,1)+time)==3
                if i==1
                    A2=kron(speye(GP(3)),kron(speye(GP(2)),Adim{i,j,1}));
                elseif i==2
                    A2=kron(speye(GP(3)),kron(Adim{i,j,1},speye(GP(1))));
                else
                    A2=kron(Adim{i,j,1},kron(speye(GP(2)),speye(GP(1))));
                end
            elseif (size(Adim,1)+time)==4
                if i==1
                    A2=kron(kron(speye(GP(4)),kron(speye(GP(3)),speye(GP(2)))),Adim{i,j,1});
                elseif i==2
                    A2=kron(kron(speye(GP(4)),kron(speye(GP(3)),Adim{i,j,1})),speye(GP(1)));
                elseif i==3
                    A2=kron(kron(speye(GP(4)),kron(Adim{i,j,1},speye(GP(2)))),speye(GP(1)));
                else
                    A2=kron(kron(Adim{i,j,1},kron(speye(GP(3)),speye(GP(2)))),speye(GP(1)));
                end
            else
                error('Not implemented yet')
            end
            A=A+Adim{i,j,2}*A2;
        end
%         'kron(kron(kron(a,b),c),d)'
    end    
end

A2=sparse(0); 
if time
    if (size(Adim,1)+time)==2
        A2=kron(At,speye(GP(1)));
    elseif (size(Adim,1)+time)==3
        A2=kron(kron(At,speye(GP(1))),speye(GP(2)));
    elseif (size(Adim,1)+time)==4
        A2=kron(kron(kron(At,speye(GP(1))),speye(GP(2))),speye(GP(3)));
    else
        error('Not implemented yet')
    end
end
if nargout==2
    if numel(A2)>1
        B=spdiags(diag(A2),0,length(A),length(A))-A;
    else
        B=-A;
    end
end
A=A2-shift1*A;


end

