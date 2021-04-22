clc;
clear;
global Nt Nx nu tau theta;
dt_max=1e-2;
theta=1/2;
tau=1.15;
nu=1e-3;
T=2;
dx=1/64;
dt=1/256;
x0=(-1:dx:1)';
%x=x0(2:end-1);
x=x0;
U0=sin(2*pi*x);
W0=U0;
Nx=length(U0);
Nx0=length(x);
Ix0=speye(Nx0);
Ix=speye(Nx);
A1=spdiags([-ones(Nx0,1),2*ones(Nx0,1),-ones(Nx0,1)]/dx^2,[-1,0,1],Nx0,Nx0);
A1(1,Nx0)=-1/dx^2;
A1(Nx0,1)=A1(1,Nx0);
A2=spdiags([-ones(Nx0,1),0*ones(Nx0,1),ones(Nx0,1)]/(2*dx),[-1,0,1],Nx0,Nx0);
A2(1,Nx0)=A2(2,1);
A2(Nx0,1)=A2(1,2);
A=nu*A1+A2;

NNt=4:2:60;
NK=length(NNt);
Err_Diag=zeros(2,NK);
CondDiag=zeros(2,NK);
TT=zeros(1,NK);
TT=zeros(1,NK);
for jn=1:NK
    Nt=NNt(jn);
    It=speye(Nt);
    dt0=dt_max/tau^(Nt-1);
    dtn=dt0*tau.^(0:Nt-1);
    t=zeros(1,Nt+1);
    for n=2:Nt+1
        t(n)=sum(dtn(1:n-1));
    end
    T=t(end);
    TT(jn)=T;
    U_ref=zeros(Nx0,Nt+1);
    U_ref(:,1)=U0;
    U_sbs=zeros(Nx0,Nt+1);
    U_sbs(:,1)=U0;
    for n=1:Nt
        U_sbs(:,n+1)=((Ix+dtn(n)*theta*A)\(Ix-dtn(n)*(1-theta)*A))*U_sbs(:,n);
        U_ref(:,n+1)=expm(-A*t(n+1))*U0; 
    end
        %% U_ParaDIAG1: solution computed by the ParaDIAG-I algorithm proposed by Martin SISC-2019
    B1=zeros(Nt,Nt);
    for n=1:Nt
        B1(n,n)=1/dtn(n);
    end
    for n=1:Nt-1
        B1(n+1,n)=-1/dtn(n+1);
    end
    B2=spdiags([(1-theta)*ones(Nt,1),theta*ones(Nt,1)],[-1,0],Nt,Nt);
    b=kron(B2\It,Ix)*[(Ix/dtn(1)-(1-theta)*A)*W0;zeros((Nt-1)*Nx,1)];
    U_ParaDIAG1=zeros(Nx,Nt+1);
    U_ParaDIAG1(:,1)=U0;
    [V,D,invV]=getV(dtn); 
    sol_stepA=StepAC(invV,b,Nx);
    sol_stepB=zeros(Nx*Nt,1);
    for n=1:Nt
         sol_stepB((n-1)*Nx+1:n*Nx,1)=(D(n)*Ix+A)\sol_stepA((n-1)*Nx+1:n*Nx,1);
    end
   sol_stepC=StepAC(V,sol_stepB,Nx);
   U_ParaDIAG1(:,2:Nt+1)=reshape(sol_stepC,Nx,Nt);
     %% U_ParaDIAG2: solution computed by the ParaDIAG-2 algorithm proposed new
    dt=T/Nt;
    [V,D0,invV,iter]=fasteigB(Nt,dt,1e-9);
    D=diag(D0);
    b=[U0/(2*dt);zeros((Nt-1)*Nx0,1)];
    U_ParaDIAG2=zeros(Nx0,Nt+1);
    U_ParaDIAG2(:,1)=U0;
    stepA=StepAC(invV,b,Nx0);
    stepB=zeros(Nx0*Nt,1);
    U_refnew=zeros(Nx0,Nt+1);
    U_refnew(:,1)=U0;
    for n=1:Nt
        stepB((n-1)*Nx0+1:n*Nx0,1)=(D(n,n)*Ix0+A)\stepA((n-1)*Nx0+1:n*Nx0,1);
        U_refnew(:,n+1)=expm(-(n*dt)*A)*U0; 
    end
   stepC=StepAC(V,stepB,Nx0);
    U_ParaDIAG2(:,2:Nt+1)=reshape(stepC,Nx0,Nt);
    
    
     Err_sbs(jn)=max(max(abs((U_ref-U_sbs))));
    Err_Diag(1,jn)=max(max(abs((U_ref-U_ParaDIAG1(1:Nx0,:)))));
    Err_Diag(2,jn)=max(max(abs((U_refnew-U_ParaDIAG2))));
    figure(2);
    semilogy(NNt(1:jn),Err_Diag(2,1:jn),'r-.s',NNt(1:jn),Err_Diag(1,1:jn),'b-.o',NNt(1:jn),Err_sbs(1:jn),'k',NNt(1:jn),TT(1:jn),'w-d','markersize',8,'linewidth',1);shg;
    set(gca,'fontname','Times New Roman','fontsize',13);
    xlabel('$N_t$','interpreter','latex','fontsize',18);
    ylabel('global error','interpreter','latex','fontsize',18);
    leg=legend('ParaDiag-I using $\Delta t$ (new)', 'ParaDiag-I using $\Delta t_n=\Delta t_{N_t}\tau^{n-N_t}$','Time-stepping TR using $\Delta t_n=\Delta t_{N_t}\tau^{n-N_t}$');
    set(leg,'interpreter','latex','fontsize',15,'location','northwest');
    %title(['$T=$',num2str(T)],'interpreter','latex','fontsize',20);
    xlim([NNt(1)-1,NNt(end)+1]);
  end
%mesh(t_ode45(1:20:length(t_ode45)),x,U_ode45(:,1:20:length(t_ode45)));shg
% nt=60;
% mesh(t(1:nt),x,U_ParaDIAG1(:,1:nt));shg
% 
% set(gca,'fontname','Times New Roman','fontsize',12);
% xlabel('time: t','fontname','Times New Roman','fontsize',20);
% ylabel('space: x','fontname','Times New Roman','fontsize',20);
% %zlabel('$\mathbf{u}_{\rm ode45}$','interpreter','latex','fontsize',20);
% xlim([0,t(nt)]);
% zlabel('$\mathbf{u}_{\rm ParaDIAG-I}$','interpreter','latex','fontsize',20);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val=StepAC(V,F,Nx)
    global Nt  
    G=zeros(Nx*Nt,1);
    for j=1:Nt
        Z=zeros(Nx,1);
        for k=1:Nt
            Z=Z+V(j,k)*F((k-1)*Nx+1:k*Nx,1);
        end
       G((j-1)*Nx+1:j*Nx,1)=Z;
    end
    val=G;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,D,invV]=getV(dtn)
    global  tau Nt;
    p=zeros(Nt-1,1);
    q=zeros(Nt-1,1);
    tD=ones(Nt,1);
    invtD=ones(Nt,1);
    for n=1:Nt-1
        p(n)=prod((1+tau.^(1:n))./(1-tau.^(1:n)));
        q(n)=(1/tau^n)*prod((1+(tau^2)*tau.^(-n:-1))./(1-tau.^(-n:-1)));
        tD(n)=1/sqrt(1+sum(abs(p(1:Nt-n)).^2));
        invtD(n)=sqrt(1+sum(abs(p(1:Nt-n)).^2));
    end
    V=toeplitz([1;p],[1,zeros(1,Nt-1)]);
    invV=toeplitz([1;q],[1,zeros(1,Nt-1)]);
    D=2./dtn;
end
