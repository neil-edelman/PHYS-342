c   = 3e8;
mu0 = pi*4e-7;

dX = 1; dY = 1; dZ = .1;
X  = 1; Y  = 0; Z  = 1;

% real space of current source
x = [X:dX:X];
y = [Y:dY:Y];
z = [-Z/2:dZ:Z/2];
% what is this? units? X:dX:X?

% current source density
Iz = 0;
Iz(1:length(x), 1, length(z)) = 0; % <- ?
Iz(1,1,:) = 1; % <- ?

% space-time for Gain calculation for radius r0
r0  = 100;
dr0 = 0.01;
dth = pi/20;
dphi= pi/20;
th  = [0:dth:pi];
phi = [0:dphi:2*pi];
w   = 1e9;
dt  = 0.1/w;
tr  = [0:dt:2*pi/w];

% iteration for the time average
PemAv = 0;
SpAv  = 0;
for itr=1:length(tr);
	Az = 0*ones(length(th),length(phi),2);
	Pem = 0;
	Bphi= 0;
	Bth = 0;

	% iteration for the surface of a sphere and 2 radii, integral for A
	for ith=1:length(th), for iphi=1:length(phi), for ir=0:1, r=ir*dr0+r0;
		% conversion to spherical coordinates
		xr = r*sin(th(ith))*cos(phi(iphi));
		yr = r*sin(th(ith))*sin(phi(iphi));
		zr = r*cos(th(ith));
		integral = 0;
		% triple for loop for volume integral
		for ix = 1:length(x), for iy = 1:length(y), for iz = 1:length(z),
			dr       = sqrt((x(ix)-xr)^2+(y(iy)-yr)^2+(z(iz)-zr)^2);
			integral = (dX*dY*dZ*Iz(ix,iy,iz)*exp(-i*w*(tr(itr)-dr/c))/(4*pi*dr))+integral;
		end; end; end;
		% in units of mu0
		Az(ith, iphi, ir + 1) = integral;
		% defining the indices for r,theta and phi in order to take derivatives
		if ir*(ith-1)*(iphi-1)>0,
			Bphi = -sin(th(ith))*(r*Az(ith,iphi,ir+1)-r0*Az(ith,iphi,ir))/(dr0*r0)-(cos(th(ith))*Az(ith,iphi,ir)-cos(th(ith-1))*Az(ith-1,iphi,ir))/(r0*dth);
			% the magnetic field component in the theta direction
			Bth  = (cos(th(ith))/(r0*sin(th(ith))))*(Az(ith,iphi,ir)-Az(ith,iphi-1,ir))/dphi;
			% Poynting vector
			S            = (c*mu0)*(real(Bphi)^2+real(Bth)^2);
			Sp(ith,iphi) = S;
			% Integration of the Poynting vector over the surface of a sphere
			Pem  = Pem+r0^2*sin(th(ith))*dth*dphi*S;
		end;
	end; end; end;

	% Time-averaged power emitted
	PemAv = PemAv+dt*w*Pem/(2*pi);
	% Calculating the time-average of the Poynting vector
	SpAv  = SpAv+dt*w*Sp/(2*pi);
end;

Pemth = mu0/(12*pi*c)*(w*Z)^2;  % Theoretical power 
Gain  = (4*pi*r0^2/PemAv)*SpAv; % Gain
Gain2 = Gain;
% fixing the first index
Gain2(1,:) = Gain(length(th),:);
Gain2(:,1) = Gain(:,length(phi));

% mesh(phi,th,Gain);xlabel('\phi [rad]');ylabel('\theta [rad]');

% 3d plot of Gain
[th0 phi0] = meshgrid(phi,th); shape=Gain2;mnz = min(min(shape));shape = shape + abs(mnz);
mns = min(min(shape));mxs = max(max(shape));rat = mxs - mns;shape=shape/rat;shape(shape<(rat/125)) = rat/125;
[x y z] = sph2cart(phi0,th0,shape);
surf(x,y,z)
hold on
colormap(flipud(jet(1024)));
set(gcf,'renderer','OpenGL','resize','off');lighting phong;shading interp
set(gcf,'color',[0 0 0]);
g = light;set(g,'Style','infinite','position',[1 0 1]);
p = light;set(p,'Style','infinite','position',[-1 0 1]);
material shiny
set(gca,'DataAspectRatio',[1 1 1])
axis equal off vis3d
box off
hold on
