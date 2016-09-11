% 2D Magnetostatic field distribution calculator using finite differences
% and Gauss-Seidel relaxation
clear  
xn=50;yn=50; %simulation size
mu=ones(xn,yn); mu(27:35,27:35)=.2; %permeability of the system
J=0*ones(xn,yn);J(20:25,20:25)=1;% Current density
% Iterative solution with periodic boundary conditions (Gauss-Seidel plus relaxation)
rel=1;
A(xn,yn)=0;A0=A; %initialization
for t=1:100,pause(0) %iteration
  for itx=1:xn, 
      for ity=1:yn, 
   A(itx,ity)=rel*(1/4)*(A(mod(itx,xn-1)+1,ity)*(1-(mu(mod(itx,xn-1)+1,ity)-mu(mod(itx-2,xn-1)+1,ity))/4/mu(itx,ity))+A(mod(itx-2,xn-1)+1,ity)*(1+(mu(mod(itx,xn-1)+1,ity)-mu(mod(itx-2,xn-1)+1,ity))/4/mu(itx,ity))+A(itx,mod(ity,yn-1)+1)*(1-(mu(itx,mod(ity,yn-1)+1)-mu(itx,mod(ity-2,yn-1)+1))/4/mu(itx,ity))+A(itx,mod(ity-2,yn-1)+1)*(1+(mu(itx,mod(ity,yn-1)+1)-mu(itx,mod(ity-2,yn-1)+1))/4/mu(itx,ity))+J(itx,ity)*mu(itx,ity))+(1-rel)*A(itx,ity);
      end;end;
 err(t)=sum(sum((A-A0).^2)); %error between iterations
 A0=A;
 end;
U3(:,:,1)=A;U3(:,:,2)=A;[Bx3,By3,Bz3]=curl(0*U3,0*U3,U3);Bx=Bx3(:,:,1);By=By3(:,:,1); %defining the curl of A
figure(1);contour(mu,'b');hold on;contour(J,'r');quiver(Bx,By);hold off %plotting the field vectors
figure(2);loglog(err);