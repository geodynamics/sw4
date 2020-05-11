%-*-octave-*--
% makeMRgrid: make a 2-D curvilinear grid with mesh refinement
%               
% function makeMRgrid( verbose )
%
function makeMRgrid( verbose )

% xMax: dimension in x
  xMax = 50;
% zMax: z-coordinate of bottom curvilinear grid
  zMax=5;
% zBot: max dimension in z
  zBot=20;
  
  if nargin < 4
    verbose=0;
  end

  Amp = 1.5;
  Ncoarse = 11; % ds0= 0.1
  xmid = 30;
  
%  Ncoarse = round(zMax/h0) + 1; % grid size of coarsest grid
  h0 = zMax/(Ncoarse-1); % corrected grid size of coarsest grid
  ds0 = 1/(Ncoarse-1);
  if (verbose)
    printf("Ncoarse=%d, h0=%e, ds0=%e\n", Ncoarse, h0, ds0);
  end

  x0 = (0:h0:xMax);
%  if (verbose) size(x0)
  k0 = 2*pi/xMax;

  topo0 = Amp*sin(k0*x0) - 0.5*Amp*cos(2*k0*x0).*exp(-(x0-xmid).^2/100);

  sref = [0.6, 0.3, 0];
  
  ds0 = 1/(Ncoarse-1);
  s0 = (sref(1):ds0:1)';
  z0 = curvigrid(s0, topo0, zMax, Ncoarse);

  h1 = 0.5*h0;
  x1 = (0:h1:xMax);

  topo1 = Amp*sin(k0*x1)- 0.5*Amp*cos(2*k0*x1).*exp(-(x1-xmid).^2/100);
  
  ds1 = 0.5*ds0;
  s1 = (sref(2):ds1:sref(1))';
  z1 = curvigrid(s1, topo1, zMax, Ncoarse);
  
  h2 = 0.5*h1;
  x2 = (0:h2:xMax);

  topo2 = Amp*sin(k0*x2)- 0.5*Amp*cos(2*k0*x2).*exp(-(x2-xmid).^2/100);
  
  ds2 = 0.5*ds1;
  s2 = (sref(3):ds2:sref(2))';
  z2 = curvigrid(s2, topo2, zMax, Ncoarse);

				% Cartesian grid
  hc = 2*h0;
  xc = (0:hc:xMax);
  Nx = size(xc,2);
  Nz = (zBot-zMax)/hc + 1;
  dsc = 1/(Nz-1);
  sc = (0:dsc:1)';
  zc = zeros(Nz,Nx);
  for k=1:Nx
    zc(:,k) = (1-sc).*zMax + sc.*(zBot);
  end

  plotgrid(xc,zc,"m");
  hold on;
  
  plotgrid(x0,z0,"b",1);
  hold on;
  plotgrid(x1,z1,"r");
  hold on;
  plotgrid(x2,z2,"k");

  axis equal
  set(gca,"fontsize",16)
  
end




