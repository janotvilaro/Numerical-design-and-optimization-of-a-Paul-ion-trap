%% sesiones 3 en adelante

%% sesion 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) Mallar los electrodos y visualizarlos
close all
clear all
delta=[0.5,2*pi/30];
R=0.01; it=7; v=[2.5,3,-1.5,1.5];
[Posmattop,Posmatbottom,Posmatring,obj]=hiperbolic(delta,R,it,v);

% Top
figure(1)
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');

% Puntos que hara de vertice
figure(2)
plot3(Posmattop(1,:),Posmattop(2,:),Posmattop(3,:),'x')
hold on
plot3(Posmatbottom(1,:),Posmattop(2,:),Posmatbottom(3,:),'x')
plot3(Posmatring(1,:),Posmatring(2,:),Posmatring(3,:),'x')

% Vertices del anillo
figure(3)
plot3(Posmatring(1,:),Posmatring(2,:),Posmatring(3,:),'x')
axis equal

% Setup final
figure(4)
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B) Hallar las cargas en los electrodos y resolver el problema 
% electrostatico
N=10;
epsilon0=1; %%%%%FIJARSE AL HACER EL FINAL
Lx=2; Ly=2; Lz=0; d=0.05;
v1 = obj.Posmat(:,obj.topol(1,:));
v2 = obj.Posmat(:,obj.topol(2,:));
v3 = obj.Posmat(:,obj.topol(3,:));
obj.cent =(v1+v2+v3)/3;
c = cross(v2-v1, v3-v1);
obj.ds = sqrt(sum(c.^2))/2;
obj.un = c./repmat(2*obj.ds,3,1); 
Z=zeros(max(size(obj.cent)));
for ii=1:max(size(obj.cent))
Z(ii,:) = int_S_1divR(obj.cent(:,ii) , v1 , v2 , v3 , obj.un , obj.cent)/(4*pi*epsilon0);
end
V0=4;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.05+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end

figure(5)
patch('Faces',obj.topol','Vertices',obj.Posmat','FaceColor','g','EdgeColor','k');

figure(6)
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
axis equal;
colormap jet; colorbar
title('Charge distribution along the meshed electrodes')

%Plot V
% [xx,yy,zz] = meshgrid(obj.cent(1,[1:30:length(obj.cent)]),obj.cent(2,[1:30:length(obj.cent)]),obj.cent(3,[1:30:length(obj.cent)]));
% coneplot(xx,yy,zz,obj.un(1,[1:30:length(obj.cent)]),obj.un(2,[1:30:length(obj.cent)]),obj.un(3,[1:30:length(obj.cent)]))%%
xs=0; ys=0; zs=0;
figure(7)
slice(xx,yy,zz, Vmat, xs,ys,zs);
colormap jet; shading interp; colorbar
axis equal;
title('Cross section of the potential produced by the hyperbolic electrodes')

%% SESIÓN 4
close all

% 1 ion

p=[0 0 0 0 1 0]; q=1; m=1; deltat=10^(-2); N=2048;   
posmat=ionprimero(p,q,m,deltat,N,Vmat,xx,yy,zz);


figure(8)
plot3(posmat(1,:),posmat(2,:),posmat(3,:))
axis([-4 4 -4 4 -4 4])
%Trayectoria sola

figure(9)
plot3(posmat(1,:),posmat(2,:),posmat(3,:))
hold on
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

% ok, vinicial=vx, y no parece desviarse on otra direccion... solo se
% acerca al electrodo.


figure(10)
comet3(posmat(1,:),posmat(2,:),posmat(3,:))
xlabel('x'); ylabel('y'); zlabel('z')


% 2 iones

%Entradas 
% x=p(1); y=p(2); z=p(3); vx=p(4); vy=p(5); vz=p(6);
% x2=p(7); y2=p(8); z2=p(9); vx2=p(10); vy2=p(11); vz2=p(12);
p=[0 0 0 1 0 0 0.1 0 0 0 0 1]; q=[-1 -1]; m=[1 1]; deltat=10^(-3); N=3000;   
[posmat1,posmat2]=dos_iones(p,q,m,deltat,N,Vmat,xx,yy,zz);

figure(11)
plot3(posmat1(1,:),posmat1(2,:),posmat1(3,:))
hold on
plot3(posmat2(1,:),posmat2(2,:),posmat2(3,:))
axis([-4 4 -4 4 -4 4])
%Trayects solas

figure(12)
plot3(posmat1(1,:),posmat1(2,:),posmat1(3,:))
hold on
plot3(posmat2(1,:),posmat2(2,:),posmat2(3,:))
hold on
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');
%jugar con num its y deltat

figure(13)
comet3(posmat1(1,:),posmat1(2,:),posmat1(3,:))
xlabel('x'); ylabel('y'); zlabel('z')
hold on
comet3(posmat2(1,:),posmat2(2,:),posmat2(3,:))


% n iones
%2 iones
%p=[x1 x2; y1 y2; z1 z2; vx1 vx2; vy1 vy2; vz1 vz2]
p=[0 0.1;0 0;0 0;1 0;0 1;0 0]; q=[-1 -1]; m=[1 1]; deltat=10^(-3); N=3500;   
posmat=n_iones(p,q,m,deltat,N,Vmat,xx,yy,zz,2);

figure(14)
plot3(posmat(1,:,1),posmat(2,:,1),posmat(3,:,1))
hold on
plot3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))
axis([-4 4 -4 4 -4 4])


figure(15)
plot3(posmat(1,:,1),posmat(2,:,1),posmat(3,:,1))
hold on
plot3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))
hold on
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

% figure(16)
% comet3(posmat(1,:,:),posmat(2,:,:),posmat(3,:,:))
% xlabel('x'); ylabel('y'); zlabel('z')
% hold on
% comet3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))


% %4 iones

%p=[x1 x2 ...; y1 y2 ...; z1 z2 ...; vx1 vx2 ...; vy1 vy2 ...; vz1 vz2 ...]
%ultima part empieza parada
p=[0 0.1 0.2 -0.1;0 0 0 0;0 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0]; q=[-1 -1 -1 -1]; 
m=[1 1 1 1]; deltat=10^(-3.5); N=3500;   
posmat=n_iones(p,q,m,deltat,N,Vmat,xx,yy,zz,4);

figure(17)
plot3(posmat(1,:,1),posmat(2,:,1),posmat(3,:,1))
hold on
plot3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))
hold on
plot3(posmat(1,:,3),posmat(2,:,3),posmat(3,:,3))
hold on
plot3(posmat(1,:,4),posmat(2,:,4),posmat(3,:,4))
axis([-4 4 -4 4 -4 4])
% trayectorías solas

%trayectorías con los electrodos
figure(18)
plot3(posmat(1,:,1),posmat(2,:,1),posmat(3,:,1))
hold on
plot3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))
hold on
plot3(posmat(1,:,3),posmat(2,:,3),posmat(3,:,3))
hold on
plot3(posmat(1,:,4),posmat(2,:,4),posmat(3,:,4))
hold on
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

%trayectorías solas pero con comet, para visualizar el movimiento
figure(19)
comet3(posmat(1,:,:),posmat(2,:,:),posmat(3,:,:))
hold on
comet3(posmat(1,:,2),posmat(2,:,2),posmat(3,:,2))
hold on
comet3(posmat(1,:,3),posmat(2,:,3),posmat(3,:,3))
hold on
comet3(posmat(1,:,4),posmat(2,:,4),posmat(3,:,4))
xlabel('x'); ylabel('y'); zlabel('z')



%% SESIÓN 5

% Pot AC con n iones. Problema adimensional sin constantes físicas

% 2 iones
global epsilon0
epsilon0=1;
p=[0 0.1;0 0;0 0;1 0;0 1;0 0]; q=[-1 -1]; m=[1 1]; deltat=10^(-4); N=2048; 
w=10^3;
posmatAC=n_ionesAC(p,q,m,deltat,N,Vmat,xx,yy,zz,2,w);

figure(20)
plot3(posmatAC(1,:,1),posmatAC(2,:,1),posmatAC(3,:,1))
hold on
plot3(posmatAC(1,:,2),posmatAC(2,:,2),posmatAC(3,:,2))
axis([-4 4 -4 4 -4 4])


figure(21)
comet3(posmatAC(1,:,:),posmatAC(2,:,:),posmatAC(3,:,:))
xlabel('x'); ylabel('y'); zlabel('z')
hold on
comet3(posmatAC(1,:,2),posmatAC(2,:,2),posmatAC(3,:,2))


% n iones

%p=[x1 x2 ...; y1 y2 ...; z1 z2 ...; vx1 vx2 ...; vy1 vy2 ...; vz1 vz2 ...]
%ultima part empieza parada
n=8;
p=rand(6,n)-rand(6,n); q=-ones(1,n);
p([1:3],:)=1.5*p([1:3],:); 
m=ones(1,n); deltat=10^(-3); N=1500; w=10^(3);
posmat=n_ionesAC(p,q,m,deltat,N,Vmat,xx,yy,zz,n,w);

for ii=1:n
figure(22)
plot3(posmat(1,:,ii),posmat(2,:,ii),posmat(3,:,ii))
hold on
%axis([-4 4 -4 4 -4 4])
end
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

%% sesión 6: ajuste de los parametros de la trampa

%% 6.1 prueba potenial muy pequeño: repulsion deberia dominar
V0=0;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.05+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end


n=2;
p=rand(6,n)-rand(6,n); q=-370*ones(1,n);
% Ponemos q muy grande para que se vea claramente la repulsion
p([1:3],:)=1.5*p([1:3],:); 
m=ones(1,n); deltat=10^(-3); N=3000; w=10^(3);
posmat=n_ionesAC(p,q,m,deltat,N,Vmat,xx,yy,zz,n,w);

for ii=1:n
figure(23)
plot3(posmat(1,:,ii),posmat(2,:,ii),posmat(3,:,ii))
hold on
%axis([-4 4 -4 4 -4 4])
end
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');



%% 6.2 prueba con frecuencia w muy baja: deberian irse los iones a los electrodos

V0=40;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.05+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end


n=2;
p=rand(6,n)-rand(6,n); q=-ones(1,n);
p([1:3],:)=1.5*p([1:3],:); 
m=ones(1,n); deltat=10^(-3); N=8000; w=0;
posmat=n_ionesAC(p,q,m,deltat,N,Vmat,xx,yy,zz,n,w);

for ii=1:n
figure(24)
plot3(posmat(1,:,ii),posmat(2,:,ii),posmat(3,:,ii))
hold on
%axis([-4 4 -4 4 -4 4])
end
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

%OBS: Parece correco, la z crece rapidamente hasta llegar a chocar contra
%el electrodo


%% 6.3 ya ponemos las constantes fisicas reales
clear all 
close all
% Vamos a diseñar una tramopa de unos 10 centimetros

delta=[0.01,2*pi/30];
R=0.002; it=7; v=[0.04,sqrt(2)*0.04,-0.03,0.03]; 
%ponemos altura 0.04 y radio sqrt(2)*altura para mas info bibliografia Paul
%trap aix marseille pag 9
[~,~,~,obj]=hiperbolic(delta,R,it,v);

figure(1)
patch('Faces',obj.topolring','Vertices',obj.Posmatring','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topoltop','Vertices',obj.Posmattop','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj.topolbottom','Vertices',obj.Posmatbottom','FaceColor','g','EdgeColor','k');

epsilon0=8.8541878176*10^(-12);
N=10;
Lx=0.05; Ly=0.05; Lz=0; d=0.001;
v1 = obj.Posmat(:,obj.topol(1,:));
v2 = obj.Posmat(:,obj.topol(2,:));
v3 = obj.Posmat(:,obj.topol(3,:));
obj.cent =(v1+v2+v3)/3;
c = cross(v2-v1, v3-v1);
obj.ds = sqrt(sum(c.^2))/2;
obj.un = c./repmat(2*obj.ds,3,1); 
Z=zeros(max(size(obj.cent)));
for ii=1:max(size(obj.cent))
Z(ii,:) = int_S_1divR(obj.cent(:,ii) , v1 , v2 , v3 , obj.un , obj.cent)/(4*pi*epsilon0);
end
V0=0.002;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end

xs=0; ys=0; zs=0;
%potencial dibujo
figure(2)
slice(xx,yy,zz, Vmat, xs,ys,zs);
colormap jet; shading interp; colorbar
axis equal;

%cargas dibujo
figure(3)
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor',[0.2 0.2 0.2]);
axis equal;
colormap jet; colorbar

% 1 ION de Berilio+ MOVIMIENTO MARIPOSA
n=1;
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 

deltat=10^(-5); 
N=50000; 
w=1000;
posmat=ionprimeroAC(p,qion,m,deltat,N,Vmat,xx,yy,zz,w);

figure(4)
plot3(posmat(1,:),posmat(2,:),posmat(3,:))
hold on
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
alpha(0.7)
axis equal;
colormap jet; colorbar

figure(5)
subplot(3,1,1)
plot([0:deltat:(N-1)*deltat],posmat(1,:))
xlabel('t (s)')
ylabel('x (m)')
subplot(3,1,2)
plot([0:deltat:(N-1)*deltat],posmat(2,:))
xlabel('t (s)')
ylabel('y (m)')
subplot(3,1,3)
plot([0:deltat:(N-1)*deltat],posmat(3,:))
xlabel('t (s)')
ylabel('z (m)')
sgtitle('Coordinates of the ion movement')
errordistancia=vecnorm(posmat);
min(errordistancia)
%% 2 iones de berilio +
% Bucle para ver que combinaciones son optimas: LARGO DE CORRER


epsilon0=8.854187817599999e-12;
% petada=zeros(length(V0),length(wvec));
petada=[];
hh=1;
for V0=[1:2:10]
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
    for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
    end
    Vmat=zeros(10,10,10);
    for jj=1:10
        for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
        end
    end

n=2;
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
wvec=[20000:10000:100000];

    for w=wvec
        deltat=(2*pi)/(50*w); 
        N=round(0.01/deltat); 
        posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
        [row, col] = find(isnan(posmat(:,:,1)));
        emptycol=isempty(col);
     if emptycol==1
         col=0;
     else
     end
        petada=[petada col(1)];
        ravg=sum(vecnorm(posmat),2);
        [ravgmax,ionmax]=max(ravg);
        str = 'V0='+string(V0)+' V'+' W='+string(w)+' rad/s';
        figure(6)
        subplot(5,9,hh)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,ionmax))); 
        xlabel('t (s)')
        ylabel('r (m)')
        title(str)
        hh=hh+1;
        end
end



%% En vista de los de arriba calculamos casos particulares

% however, if the charge of the ions is
% increased (or the trap potential is reduced), their mutual
% interaction should have a greater impact.

V0=3;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*8.854187817599999e-12)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end


n=2;
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.02*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
deltat=1*10^(-6); 
N=100000; 
w=50000;
posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
for ii=1:n
figure(7)
plot3(posmat(1,1:end-1,ii),posmat(2,1:end-1,ii),posmat(3,1:end-1,ii))
hold on
% plot3(posmat(1,end,ii),posmat(2,end,ii),posmat(3,end,ii),'x','MarkerSize',200)
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
alpha(0.2)
axis equal;
colormap jet; colorbar
end

figure(8)
subplot(3,1,1)
plot([0:deltat:(N-1)*deltat],posmat(1,:,1)); hold on
plot([0:deltat:(N-1)*deltat],posmat(1,:,2));
xlabel('t (s)')
ylabel('x (m)')
subplot(3,1,2)
plot([0:deltat:(N-1)*deltat],posmat(2,:,1)); hold on
plot([0:deltat:(N-1)*deltat],posmat(2,:,2));
xlabel('t (s)')
ylabel('y (m)')
subplot(3,1,3)
plot([0:deltat:(N-1)*deltat],posmat(3,:,1)); hold on
plot([0:deltat:(N-1)*deltat],posmat(3,:,2))
xlabel('t (s)')
ylabel('z (m)')
sgtitle('Coordinates of the ion movement')

%% Optimizacion para N iones
epsilon0=8.854187817599999e-12;
% petada=zeros(length(V0),length(wvec));
petada=[];
hh=1;
for n=[20]
for V0=[3 5 7]
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
    for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
    end
    Vmat=zeros(10,10,10);
    for jj=1:10
        for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
        end
    end

p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
wvec=[50000 62500 75000];

    for w=wvec
        deltat=(2*pi)/(50*w); 
        N=round(0.01/deltat); 
        posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
        [row, col] = find(isnan(posmat(:,:,1)));
        emptycol=isempty(col);
     if emptycol==1
         col=0;
     else
     end
        petada=[petada col(1)];
        ravg=sum(vecnorm(posmat),2);
        [ravgmax,ionmax]=max(ravg);
        str = 'N='+string(n)+' V0='+string(V0)+' V'+' W='+string(w)+' rad/s';
        figure(15)
        subplot(4,3,hh)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,ionmax))); 
        xlabel('t (s)')
        ylabel('r (m)')
        title(str)
        hh=hh+1;
        end
end
end
%% N iones
% LARGO DE CORRER (5 mins max)

V0=7;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*8.854187817599999e-12)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end


kk=1;
for n=[5 7 10 15 20]
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.02*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
deltat=1*10^(-6); 
N=10000; 
w=62500;
posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
for ii=1:n
figure(9+kk)
plot3(posmat(1,1:end-1,ii),posmat(2,1:end-1,ii),posmat(3,1:end-1,ii))
hold on
% plot3(posmat(1,end,ii),posmat(2,end,ii),posmat(3,end,ii),'x','MarkerSize',200)
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
alpha(0.2)
axis equal;
colormap jet; colorbar
end
ravg=sum(vecnorm(posmat),2);
[ravgmax,ionmax]=max(ravg);
str = 'N='+string(n);
figure(9)
subplot(3,2,kk)
plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,ionmax))); 
xlabel('t (s)')
ylabel('r (m)')
title(str)
kk=kk+1;
end

n=20;
w=62500;
V0=7;
N=100000; 
posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
for ii=1:n
figure(16)
plot3(posmat(1,1:end-1,ii),posmat(2,1:end-1,ii),posmat(3,1:end-1,ii))
hold on
% plot3(posmat(1,end,ii),posmat(2,end,ii),posmat(3,end,ii),'x','MarkerSize',200)
patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
alpha(0.2)
axis equal;
colormap jet; colorbar
end
ravg=sum(vecnorm(posmat),2);
[ravgmax,ionmax]=max(ravg);
str = 'N='+string(n)+' V0='+string(V0)+' V'+' W='+string(w)+' rad/s';
figure(17)
plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,ionmax))); 
xlabel('t (s)')
ylabel('r (m)')
title(str)

%% Conclusiones Extra
%% Vamos a confinar 1 ion al maximo: quantum computing.
n=1;
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
pppp=1;
minerrordistancia=[];
for V0=[0.001 0.01 0.1 1 10]
    b=V0/2*ones(max(size(obj.cent)),1);
    b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
    Nfb1p=numel(obj.ds);
    q=Z\b;
    % Computacion de V y E
    x = linspace(-Lx,Lx,10);
    y = linspace(-Ly,Ly,10);
    z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
    [xx,yy,zz] = meshgrid(x,y,z);
    V=zeros(10*10*10,1);
    for ii=1:max(size(obj.cent))
        RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los
        %puntos de carga y los del meshgrid.
        RR=vecnorm(RR);
        V=V+1/(4*pi*8.854187817599999e-12)*q(ii)*obj.ds(ii)./RR'; %Coulomb
    end
    Vmat=zeros(10,10,10);
    for jj=1:10
        for kk=1:10
            Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
        end
    end
    for w=[100 1000 10000 50000]
        deltat=(2*pi)/(50*w);
        N=round(0.1/deltat);
        posmat=ionprimeroAC(p,qion,m,deltat,N,Vmat,xx,yy,zz,w);
        figure(18)
        subplot(5,4,pppp)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat))
        xlabel('t (s)')
        ylabel('r (m)')
        sgtitle('1 ion confinement for different settings')
        str = 'V0='+string(V0)+' V'+' W='+string(w)+' rad/s';
        title(str)
        errordistancia=vecnorm(posmat);
        minerrordistancia=[minerrordistancia min(errordistancia)];
        pppp=pppp+1;
        figure(20+pppp)
        plot3(posmat(1,:),posmat(2,:),posmat(3,:))
        hold on
        patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor','none');
        alpha(0.7)
        axis equal;
        colormap jet; colorbar
    end
end



%% Jugar con masas y cargas relativas con dos iones.
% Con el setup optimo hallado para mrel=qrel=1, vemos el efecto de
% augmenatr y disminuir qrel o mrel.

V0=3;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*8.854187817599999e-12)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end

n=2;
p=rand(6,n)-rand(6,n); 
p([1:3],:)=0.02*p([1:3],:); 
deltat=1*10^(-6); 
N=10000; 
w=50000;
ppp=1;
for mrel=2*[0.25 0.5 0.75]
    for qrel=[2]            
        m=9.012182*1.6605*10^(-27)*ones(1,n); 
        qion=1.602*10^(-19)*ones(1,n);
        m=[m(1) mrel*m(2)];
        qion=[qion(1) qrel*qion(2)];
        posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
        ravg=sum(vecnorm(posmat),2);
        [ravgmax,ionmax]=max(ravg);
        figure(19)
        subplot(3,3,ppp)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,1)))
        hold on
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,2)))
        xlabel('t (s)')
        ylabel('r (m)')
        sgtitle('2 ion confinement for different settings')
        str = 'qrel='+string(qrel)+' mrel='+string(mrel);
        title(str)
%             errordistancia=vecnorm(posmat);
%             minerrordistancia=[minerrordistancia min(errordistancia)];
        ppp=ppp+1;
    end
end

for mrel=1000*[0.0025 0.005 0.0075]
    for qrel=[10]            
        m=9.012182*1.6605*10^(-27)*ones(1,n); 
        qion=1.602*10^(-19)*ones(1,n);
        m=[m(1) mrel*m(2)];
        qion=[qion(1) qrel*qion(2)];
        posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
        ravg=sum(vecnorm(posmat),2);
        [ravgmax,ionmax]=max(ravg);
        figure(19)
        subplot(3,3,ppp)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,1)))
        hold on
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,2)))
        xlabel('t (s)')
        ylabel('r (m)')
        str = 'qrel='+string(qrel)+' mrel='+string(mrel);
        title(str)
%             errordistancia=vecnorm(posmat);
%             minerrordistancia=[minerrordistancia min(errordistancia)];
        ppp=ppp+1;
    end
end

for mrel=10000*[0.0025 0.005 0.0075]
    for qrel=[100]            
        m=9.012182*1.6605*10^(-27)*ones(1,n); 
        qion=1.602*10^(-19)*ones(1,n);
        m=[m(1) mrel*m(2)];
        qion=[qion(1) qrel*qion(2)];
        posmat=n_ionesAC_SI(p,qion,m,deltat,N,Vmat,xx,yy,zz,n,w);
        ravg=sum(vecnorm(posmat),2);
        [ravgmax,ionmax]=max(ravg);
        figure(19)
        subplot(3,3,ppp)
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,1)))
        hold on
        plot([0:deltat:(N-1)*deltat],vecnorm(posmat(:,:,2)))
        xlabel('t (s)')
        ylabel('r (m)')
        str = 'qrel='+string(qrel)+' mrel='+string(mrel);
        title(str)
%             errordistancia=vecnorm(posmat);
%             minerrordistancia=[minerrordistancia min(errordistancia)];
        ppp=ppp+1;
    end
end