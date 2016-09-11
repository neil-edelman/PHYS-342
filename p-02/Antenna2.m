c=3e8;mu0=pi*4e-7; %constants
dX=1;dY=1;dZ=.001;X=0;Y=0;Z=0.1;x=[X:dX:X];y=[Y:dY:Y];z=[-Z/2:dZ:Z/2]; % real space of current source
Iz=0;Iz(1,1,length(z))=0;Iz(1,1,:)=ones(length(z),1); % current source density
r0=100;dr0=0.01;dth=0.2;dphi=0.2;th=[0:dth:pi];phi=[0:dphi:2*pi];w=1e9;dt=0.1/w;tr=[0:dt:2*pi/w];  % space-time for Gain calculation for radius r0
PemAv=0;SpAv=0;for itr=1:length(tr); % iteration for the time average
Az=0*ones(length(th),length(phi),2);Pem=0;Bphi=0;Bth=0;
for ith=1:length(th),for iphi=1:length(phi),for ir=0:1, r=ir*dr0+r0; % iteration for the surface of a sphere and 2 radii (
xr=r*sin(th(ith))*cos(phi(iphi));yr=r*sin(th(ith))*sin(phi(iphi));zr=r*cos(th(ith)); % conversion to spherical coordinates
                % integral for A
integral=0;for ix=1:length(x),for iy=1:length(y),for iz=1:length(z), % triple for loop for volume integral
dr=sqrt((x(ix)-xr)^2+(y(iy)-yr)^2+(z(iz)-zr)^2);
integral=(dX*dY*dZ*Iz(ix,iy,iz)*exp(-i*w*(tr(itr)-dr/c))/(4*pi*dr))+integral;
        end;end;end;
Az(ith,iphi,ir+1)=integral; % in units of mu0
if ir*(ith-1)*(iphi-1)>0, % defining the indices for r,theta and phi in order to take derivatives
Bphi=-sin(th(ith))*(r*Az(ith,iphi,ir+1)-r0*Az(ith,iphi,ir))/(dr0*r0)-(cos(th(ith))*Az(ith,iphi,ir)-cos(th(ith-1))*Az(ith-1,iphi,ir))/(r0*dth);
Bth=(cos(th(ith))/(r0*sin(th(ith))))*(Az(ith,iphi,ir)-Az(ith,iphi-1,ir))/dphi; % the magnetic field component in the theta direction
S=(c*mu0)*(real(Bphi)^2+real(Bth)^2);Sp(ith,iphi)=S; % Poynting vector
Pem=Pem+r0^2*sin(th(ith))*dth*dphi*S; % Integration of the Poynting vector over the surface of a sphere
end;
end;end;end;
PemAv=PemAv+dt*w*Pem/(2*pi); % Time-averaged power emitted
SpAv=SpAv+dt*w*Sp/(2*pi); % Calculating the time-average of the Poynting vector
end;
Pemth=mu0/(12*pi*c)*(w*Z)^2; % Theoretical power 
Gain=(4*pi*r0^2/PemAv)*SpAv; % Gain
mesh(phi,th,Gain);xlabel('\phi [rad]');ylabel('\theta [rad]');