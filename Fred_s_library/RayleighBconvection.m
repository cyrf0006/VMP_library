%Simulating Rayleigh-Benard convection in thin layer of fluid(2D)
...Bottom and top walls are rigid, impermeable and isothermal with no-slip
...Side walls are adiabatic (insulated)
...Boussinesq approximation is used to model the buoyancy driven flow
%For convection rolls, one set of conditions are:
...H=0.2,L=0.82,TN=25,TS=70,Re=1e2,Pr=7,Ra=1e2

%%
%Specifying Parameters
clear all
close all
nx=100;                              %Number of steps in space(x)
ny=70;                               %Number of steps in space(y)
tf=1e3;                              %Final time
dt=1e-2;                             %Width of time step
nt=ceil(tf/dt);                      %Number of time steps
dt=tf/nt;                            %Corrected time step
H=2;                                 %Height of the container
L=5;                                 %Length of the container
dx=L/(nx-1);                         %Width of space step(x)
dy=H/(ny-1);                         %Width of space step(y)
x=0:dx:L;                            %Range of x(0,2) and specifying grid points
y=0:dy:H;                            %Range of y(0,5) and specifying grid points
TN=10;                               %Top wall temperature
TS=9;                               %Bottom wall temperature
u=zeros(nx+1,ny);                    %Preallocating u
v=zeros(nx,ny+1);                    %Preallocating v
p=zeros(nx,ny);                      %Preallocating p
S=zeros(nx,ny);                      %Preallocating S
uplot=zeros(nx,ny);                  %Preallocating uplot
vplot=zeros(nx,ny);                  %Preallocating vplot
To=min(TS,TN);                       %Initial temperature
T=To*ones(nx,ny)*(ones(70,1).*(1:70)*1\70);                    %Preallocating T (Initial conditions)
T(:, 1:35) = TN;
Re=1e2;                              %Reynolds number
Pr=7;                                %Prandtl number
Pe=Re*Pr;                            %Peclet number
Ra=1e2;                              %Rayleigh number
Gr=Ra/Pr;                            %Grashoff number
Tstar=T; 
ustar=u; uhalf=u; uconv=u;
vstar=v; vhalf=v; vconv=v;
TnE=0;   TnW=0;
UN=0;    VN=0;      
US=0;    VS=0;
UE=0;    VW=0;                  
UW=0;    VE=0;



%keyboard


%%
%IMPLICIT DIFFUSION:-
%B.C vector(for u)
bcu=zeros(nx-1,ny-2);
bcu(1,:)=2*UW/dx^2; bcu(end,:)=2*UE/dx^2;  %Dirichlet B.Cs
bcu(:,1)=US/dy^2; bcu(:,end)=UN/dy^2;      %Dirichlet B.Cs
%B.Cs at the corners:
bcu(1,1)=2*UW/dx^2+US/dy^2; bcu(end,1)=2*UE/dx^2+US/dy^2;
bcu(1,end)=2*UW/dx^2+UN/dy^2; bcu(end,end)=2*UE/dx^2+UN/dy^2;
bcu=dt*bcu/Re;

%B.C vector(for v)
bcv=zeros(nx-2,ny-1);
bcv(1,:)=VW/dx^2; bcv(end,:)=VE/dx^2;      %Dirichlet B.Cs
bcv(:,1)=2*VS/dy^2; bcv(:,end)=2*VN/dy^2;  %Dirichlet B.Cs
%B.Cs at the corners:
bcv(1,1)=VW/dx^2+2*VS/dy^2; bcv(end,1)=VE/dx^2+2*VS/dy^2;
bcv(1,end)=VW/dx^2+2*VN/dy^2; bcv(end,end)=VE/dx^2+2*VN/dy^2;
bcv=dt*bcv/Re;

%Central difference operator(for u)
e=ones(nx-1,1);i=ones(ny-2,1);
Ax=spdiags(e*[1 -2 1],-1:1,nx-1,nx-1);      
Ay=spdiags(i*[1 -2 1],-1:1,ny-2,ny-2);      
Ax(1,1)=-3;Ax(end,end)=-3;
A=kron(Ay/dy^2,speye(nx-1))+kron(speye(ny-2),Ax/dx^2);
Du=speye((nx-1)*(ny-2))-dt*A/Re;
pu=symamd(Du);[Lu Uu]=lu(Du(pu,pu));

%Central difference operator(for v)
e=ones(nx-2,1);i=ones(ny-1,1);
Ax=spdiags(e*[1 -2 1],-1:1,nx-2,nx-2);      
Ay=spdiags(i*[1 -2 1],-1:1,ny-1,ny-1);      
Ay(1,1)=-3;Ay(end,end)=-3;
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-1),Ax/dx^2);
Dv=speye((nx-2)*(ny-1))-dt*A/Re;
pv=symamd(Dv);[Lv Uv]=lu(Dv(pv,pv));

%%
%Calculating the coefficient matrix for the PPE
e=ones(nx-2,1);i=ones(ny-2,1);
Ax=spdiags(e*[1 -2 1],-1:1,nx-2,nx-2);      
Ay=spdiags(i*[1 -2 1],-1:1,ny-2,ny-2); 
Ax(1,1)=-1;Ax(end,end)=-1;
Ay(1,1)=-1;Ay(end,end)=-1;
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
pp=symamd(A);[Lp Up]=lu(A(pp,pp));

%%
%B.C vector(for T)
bcT=zeros(nx-2,ny-2);
bcT(1,:)=-TnW/dx; bcT(end,:)=TnE/dx;       %Neumann B.Cs
bcT(:,1)=TS/dy^2; bcT(:,end)=TN/dy^2;      %Dirichlet B.Cs
%B.Cs at the corners:
bcT(1,1)=-TnW/dx+TS/dy^2; bcT(end,1)=TnE/dx+TS/dy^2;
bcT(1,end)=-TnW/dx+TN/dy^2; bcT(end,end)=TnE/dx+TN/dy^2;
bcT=dt*bcT/Pe;

%Central difference operator(for T)
e=ones(nx-2,1);i=ones(ny-2,1);
Tx=spdiags(e*[1 -2 1],-1:1,nx-2,nx-2);      
Ty=spdiags(i*[1 -2 1],-1:1,ny-2,ny-2);
Tx(1,1)=-1; Tx(end,end)=-1;
Tt=kron(Ty/dy^2,speye(nx-2))+kron(speye(ny-2),Tx/dx^2);
Tt=speye((nx-2)*(ny-2))-dt*Tt/Pe;
pt=symamd(Tt);[Lt Ut]=lu(Tt(pt,pt));

%%
%Boundary conditions
T(1,:)=T(2,:)-TnW*dx; T(end,:)=T(end-1,:)+TnE*dx;T(:,1)=TS; T(:,end)=TN;
u(1,:)=2*UW-u(2,:); u(end,:)=2*UE-u(end-1,:); u(:,1)=US; u(:,end)=UN;
v(1,:)=VW; v(end,:)=VE; v(:,1)=2*VS-v(:,2); v(:,end)=2*VN-v(:,end-1);
keyboard
%%
%Evaluating temperature and velocity field at each time step
i=2:nx-1;
j=2:ny-1;
for it=0:nt
    uplot(1:nx,1:ny)=0.5*(u(1:nx,1:ny)+u(2:nx+1,1:ny));
    vplot(1:nx,1:ny)=0.5*(v(1:nx,1:ny)+v(1:nx,2:ny+1));
    
    if(rem(it,30)==0)        
        quiver(x,y,uplot',vplot',.6,'k');
        axis equal
        axis([0 L 0 H])
        hold on
        pcolor(x,y,T');
        colormap(jet)
        colorbar
        shading interp
        axis equal
        axis([0 L 0 H])
        title({['Rayleigh-Benard Convection with Pe = ',num2str(Pe),' and Gr = ',num2str(Gr)];['time(\itt) = ',num2str(dt*it)]})
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('Spatial co-ordinate (y) \rightarrow')
        drawnow;
        hold off
    end
    
    Tn=T;
    Tstar(i,j)=Tn(i,j)-(dt/dx/2)*(u(i+1,j).*(Tn(i,j)+Tn(i+1,j))-u(i,j).*(Tn(i,j)+Tn(i-1,j)))...
        -(dt/dy/2)*(v(i,j+1).*(Tn(i,j)+Tn(i,j+1))-v(i,j).*(Tn(i,j)+Tn(i,j-1)));
    %Thermal diffusion:
    %(1) Explicit central difference
    %{
    T(i,j)=Tstar(i,j)+(dt/Pe/dx^2)*(Tn(i+1,j)-2*Tn(i,j)+Tn(i-1,j))...
        +(dt/Pe/dy^2)*(Tn(i,j+1)-2*Tn(i,j)+Tn(i,j-1));
    %}
    %(2) Implicit central difference
    
    t=reshape(Tstar(2:end-1,2:end-1)+bcT,[],1);
    t(pt)=Ut\(Lt\t(pt));
    t=reshape(t,nx-2,ny-2);
    T(2:end-1,2:end-1)=t;
    %}
    T(1,:)=T(2,:)-TnW*dx; T(end,:)=T(end-1,:)+TnE*dx;
    T(:,1)=TS; T(:,end)=TN;
    
    un=u; vn=v;
    %Convective terms(hyperbolic):
    %(1) Central differencing with explicit Euler time march
    %{
    uconv(i+1,j)=un(i+1,j)-(dt/(4*dx))*((un(i+2,j)+un(i+1,j)).^2-(un(i+1,j)+un(i,j)).^2)...
        -(dt/(4*dy))*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-(un(i+1,j)+un(i+1,j-1)).*(vn(i,j)+vn(i+1,j)));
    vconv(i,j+1)=vn(i,j+1)-(dt/(4*dx))*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-(un(i,j)+un(i,j-1)).*(vn(i,j+1)+vn(i-1,j+1)))...
        -(dt/(4*dy))*((vn(i,j+2)+vn(i,j+1)).^2-(vn(i,j+1)+vn(i,j)).^2);
    %}
    %(2) Mac-Cormack method
    %{
    %Predictor step(Forward difference)
    uhalf(i+1,j)=un(i+1,j)-(dt/dx/2)*((un(i+1,j)+un(i+2,j)).^2-4*un(i+1,j).^2)...
        -(dt/dy/2)*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-un(i+1,j).*(vn(i,j)+vn(i,j+1)+vn(i+1,j+1)+vn(i+1,j)));
    vhalf(i,j+1)=vn(i,j+1)-(dt/dx/2)*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-vn(i,j+1).*(un(i,j+1)+un(i,j)+un(i+1,j)+un(i+1,j+1)))...
        -(dt/dy/2)*((vn(i,j+1)+vn(i,j+2)).^2-4*vn(i,j+1).^2);
    %Corrector step(Backward difference)
    uconv(i+1,j)=0.5*(un(i+1,j)+uhalf(i+1,j))-(dt/dx/4)*(4*uhalf(i+1,j).^2-(uhalf(i+1,j)+uhalf(i,j)).^2)...
        -(dt/dy/4)*(uhalf(i+1,j).*(vhalf(i,j)+vhalf(i,j+1)+vhalf(i+1,j+1)+vhalf(i+1,j))-(uhalf(i+1,j)+uhalf(i+1,j-1)).*(vhalf(i,j)+vhalf(i+1,j)));
    vconv(i,j+1)=0.5*(vn(i,j+1)+vhalf(i,j+1))-(dt/dx/4)*(vhalf(i,j+1).*(uhalf(i,j+1)+uhalf(i,j)+uhalf(i+1,j)+uhalf(i+1,j+1))-(uhalf(i,j)+uhalf(i,j+1)).*(vhalf(i,j+1)+vhalf(i-1,j+1)))...
        -(dt/dy/4)*(4*vhalf(i,j+1).^2-(vhalf(i,j)+vhalf(i,j+1)).^2);
    %}
    %(3) Richtmyer method
    
    %Predictor step(Lax-Friedrich)
    uhalf(i+1,j)=0.5*(un(i+2,j)+un(i,j))-(dt/dx/8)*((un(i+2,j)+un(i+1,j)).^2-(un(i+1,j)+un(i,j)).^2)...
        -(dt/dy/8)*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-(un(i+1,j)+un(i+1,j-1)).*(vn(i,j)+vn(i+1,j)));
    vhalf(i,j+1)=0.5*(vn(i,j+2)+vn(i,j))-(dt/dx/8)*((un(i+1,j)+un(i+1,j+1)).*(vn(i,j+1)+vn(i+1,j+1))-(un(i,j)+un(i,j-1)).*(vn(i,j+1)+vn(i-1,j+1)))...
        -(dt/dy/8)*((vn(i,j+2)+vn(i,j+1)).^2-(vn(i,j+1)+vn(i,j)).^2);
    %Corrector step(Leapfrog)
    uconv(i+1,j)=un(i+1,j)-(dt/(4*dx))*((uhalf(i+2,j)+uhalf(i+1,j)).^2-(uhalf(i+1,j)+uhalf(i,j)).^2)...
        -(dt/(4*dy))*((uhalf(i+1,j)+uhalf(i+1,j+1)).*(vhalf(i,j+1)+vhalf(i+1,j+1))-(uhalf(i+1,j)+uhalf(i+1,j-1)).*(vhalf(i,j)+vhalf(i+1,j)));
    vconv(i,j+1)=vn(i,j+1)-(dt/(4*dx))*((uhalf(i+1,j)+uhalf(i+1,j+1)).*(vhalf(i,j+1)+vhalf(i+1,j+1))-(uhalf(i,j)+uhalf(i,j-1)).*(vhalf(i,j+1)+vhalf(i-1,j+1)))...
        -(dt/(4*dy))*((vhalf(i,j+2)+vhalf(i,j+1)).^2-(vhalf(i,j+1)+vhalf(i,j)).^2);
    %}
    %Buoyancy term (Boussinesq aproximation):
    vconv(i,j+1)=vconv(i,j+1)+(Gr/Re^2)*(0.5*(T(i,j)+T(i,j+1))-To);
    
    %Diffusive terms(parabolic):
    %(1) Explicit central difference (spatial)
    %{
    ustar(i+1,j)=uconv(i+1,j)+(dt/Re/dx^2)*(un(i+2,j)-2*un(i+1,j)+un(i,j))...
        +(dt/Re/dy^2)*(un(i+1,j+1)-2*un(i+1,j)+un(i+1,j-1));
    vstar(i,j+1)=vconv(i,j+1)+(dt/Re/dx^2)*(vn(i+1,j+1)-2*vn(i,j+1)+vn(i-1,j+1))...
        +(dt/Re/dy^2)*(vn(i,j+2)-2*vn(i,j+1)+vn(i,j));
    %}
    %(2) Implicit central difference (spatial)
    
    U=reshape(uconv(2:end-1,2:end-1)+bcu,[],1);
    U(pu)=Uu\(Lu\U(pu));
    U=reshape(U,nx-1,ny-2);
    ustar(2:end-1,2:end-1)=U;
    V=reshape(vconv(2:end-1,2:end-1)+bcv,[],1);
    V(pv)=Uv\(Lv\V(pv));
    V=reshape(V,nx-2,ny-1);
    vstar(2:end-1,2:end-1)=V;
    %}     
    %Pressure Poisson equation(elliptic):
    S(i,j)=(ustar(i+1,j)-ustar(i,j))/dx...
        +(vstar(i,j+1)-vstar(i,j))/dy;
    s=reshape(S(2:end-1,2:end-1),[],1);
    s(pp)=Up\(Lp\s(pp));
    s=reshape(s,nx-2,ny-2);
    p(2:end-1,2:end-1)=s;
    p(1,:)=p(2,:);p(end,:)=p(end-1,:);
    p(:,1)=p(:,2);p(:,end)=p(:,end-1);
    
    %Velocity field update at time t+dt, along with pressure correction
    u(i+1,j)=ustar(i+1,j)-(p(i+1,j)-p(i,j))/dx;
    v(i,j+1)=vstar(i,j+1)-(p(i,j+1)-p(i,j))/dy;
    u(1,:)=2*UW-u(2,:); u(end,:)=2*UE-u(end-1,:); u(:,1)=US; u(:,end)=UN;
    v(1,:)=VW; v(end,:)=VE; v(:,1)=2*VS-v(:,2); v(:,end)=2*VN-v(:,end-1); ...
           
    keyboard
end
