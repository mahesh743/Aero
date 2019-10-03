% vortex panel method for four digit NACA series in general, using Warren F. Philip
% Mahesh Surendran
clc
clear
close all
format long


% initialization 
n = 100; % number of nodes
c = 1; % chord lenght
V_inf =1; % free stream in m/s
AOA =-2.1; % angle of attack 
alpha= AOA*pi/180;
dtheta = 2*pi/(n-1); % delta theta 

%NACA(mpxx)4 digit series
xx = 12/100; %  precentage of thckness 
m  = 0/100; % precentage maximum camberline
p =  0/100;  % location of maximum thinkness

% cosine clustering 
for i = 1:n/2
    Xn(n/2 + i)= c*(0.5)*(1-cos((i-0.5)*dtheta));
    Xn(n/2 + 1 -i)=c*(0.5)*(1-cos((i-0.5)*dtheta));
    
end

% calculating half thickness at a given x co-ordinate
for i = 1:n/2
yt(i)= 5*xx*(0.2969*sqrt(Xn(i))-0.1260*Xn(i)-0.3516*Xn(i)^2+0.2843*Xn(i)^3-0.1036*Xn(i)^4);

end


% calculating mean camber length yc and adding half-thickness yt to it. 
for i =1:n/2
    if (Xn(i)>=0 && Xn(i)<=p*c)
        yc(i) = c*m/p^2*(2*p*(Xn(i)/c)-(Xn(i)/c)^2);
    else
        yc(i) = c*m/(1-p)^2*((1-2*p)+2*p*(Xn(i)/c)-(Xn(i)/c)^2);     
    end
end

% calculation of thetas 
for i = 1:n/2
     if (Xn(i)>=0 && Xn(i)<=p*c)
         dycdx(i)=2*m/p^2*(p-(Xn(i)/c));
     else 
         dycdx(i)=2*m/(1-p)^2*(p-(Xn(i)/c));
     end
     theta(i) = atan(dycdx(i));
end


% Nodal points i.e X and Y
for i = 1:n/2
    Yn(i)= yc(i)-yt(i)*(cos(theta(i))); % lower
    Yn(n+1-i)= yc(i)+yt(i)*(cos(theta(i)));% upper
    
    Xn(i)= Xn(i) -(yt(i)*sin(theta(i)));
    Xn(n+1-i)= Xn(n+1-i)+(yt(i)*sin(theta(i)));
end

% panel centers
for i = 1:n-1
    Xc(i)= (Xn(i)+Xn(i+1))/2;
    Yc(i)= (Yn(i)+Yn(i+1))/2;
end
A=zeros(n,n);
 % airfoil coordinates
for i = 1:n-1 % looping through the panels there are n-1 panels
    x=Xc(i);
    y=Yc(i);
    
    for j=1:n-1 % influence of all other panels on the ith panel 
        
        % Airfoil coordinates
        l(j)=sqrt(((Xn(j+1))-Xn(j))^2+(Yn(j+1)-Yn(j))^2);
        
        eta=1/l(j)*(((Xn(j+1))-Xn(j))*(x-Xn(j)) + ((Yn(j+1)-Yn(j))*(y-Yn(j))));
        neta=1/l(j)*(-(Yn(j+1)-Yn(j))*(x-Xn(j)) + ((Xn(j+1))-Xn(j))*(y-Yn(j)));
        phi= atan2(neta*l(j),neta^2+eta^2-eta*l(j));
        xi=0.5*log((eta^2+neta^2)/((eta-l(j))^2+neta^2));
        
        %Panel Coefficient Matrix
        P11=(1/(2*pi*l(j)^2))*((Xn(j+1)-Xn(j))*((l(j)-eta)*phi+neta*xi)+ (-(Yn(j+1)-Yn(j)))*((neta*phi-(l(j)-eta)*xi-l(j))));
        P12=(1/(2*pi*l(j)^2))*((Xn(j+1)-Xn(j))*(eta*phi-neta*xi) + (-(Yn(j+1)-Yn(j)))*(-neta*phi-eta*xi+l(j))); 
        P21=(1/(2*pi*l(j)^2))*((Yn(j+1)-Yn(j))*((l(j)-eta)*phi+neta*xi) + ((Xn(j+1)-Xn(j))*((neta*phi-(l(j)-eta)*xi-l(j)))));
        P22=(1/(2*pi*l(j)^2))*((Yn(j+1)-Yn(j))*(eta*phi-neta*xi) + ((Xn(j+1)-Xn(j))*(-neta*phi-eta*xi+l(j))));
        
        % Airfoil coeffcient matrix
        A(i,j)=   A(i,j)+(((Xn(i+1)-Xn(i))/l(i))*P21 -((Yn(i+1)-Yn(i))/l(i))*P11);
        A(i,j+1)= A(i,j+1)+(((Xn(i+1)-Xn(i))/l(i))*P22 -((Yn(i+1)-Yn(i))/l(i))*P12);
        
    end
    
end
 A(n,1)=1;
 A(n,n)=1;
 
for i=1:n-1

    bla(i)=(((Yn(i+1)-Yn(i))*cos(alpha))-((Xn(i+1)-Xn(i))*sin(alpha)))/l(i);
end
    bla(n)=0;


% inverse multiplication             
gamma=inv(A)*(V_inf.*bla');
%gamma=A/(V_inf.*bla);

% Recovering velocities at surfaces
for i = 1:n-1
    V(i)=(gamma(i,1)+gamma(i+1,1))/2;
end


Cp= 1-(V.^2/V_inf^2);

Cl=0;

for i =1:n-1
    cl(i) =l(i)*((gamma(i)+gamma(i+1))/V_inf);
    Cmle(i)=l(i)*((2*Xn(i)*gamma(i)+Xn(i)*gamma(i+1)+Xn(i+1)*gamma(i)+2*Xn(i+1)*gamma(i+1))*cos(alpha)/V_inf...
           +(2*Yn(i)*gamma(i)+Yn(i)*gamma(i+1)+Yn(i+1)*gamma(i)+2*Yn(i+1)*gamma(i+1))*sin(alpha)/V_inf);
       
    
end
Cmle=-1/3*sum(Cmle)
cl=sum(cl)

% symmetric airfoil solution, NACA 0012 can be considered thin airfoil 
%  Ccmle=-pi/2*alpha
%  Ccl=2*pi*alpha



% Post processing
plot(Xn,Yn)
hold on
plot(Xn(1:n/2),yc)
plot(Xn(1:n/2),yt)
plot(Xc,Yc,'*')
title('Airfoil NACA 2412')
legend('2412','Camber line','Thickness','Control points')
axis equal

% Cp around the airfoil
figure(2)
plot(Xc,Cp)
hold on
set(gca,'Ydir','reverse')
title('Cp Vs X-distance')
xlabel('X-distance')
ylabel('Cp')
legend('Coefficient of pressure')
