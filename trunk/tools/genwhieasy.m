%
% Generate whi-file from Geodyn data.
%
% genwhieasy( filename, geoname, geodir, geon, geom, freq, ...
%		         dt, L, D, hbox )
%
%        Input: filename - Name of whi data file to create.
%               geoname  - Names of Geodyn station files are
%                           geoname_i_j
%               geodir   - Directory of Geodyn station files.
%               geon, geom - Size of Geodyn station files grid.
%               freq     - Filter data in time to this frequency.
%               dt       - Output data with this uniform time step.
%               L, D     - Output data on sides of a box of 
%                          size [-L/2,L/2]x[-L/2,L/2]x[0,D]
%               hbox     - Use this resolution (grid spacing) for the output data.
%
function genwhieasy( filename, geoname, geodir, geon, geom, freq, ...
		     dt, L, D, hbox )

n = length(geodir);
if geodir(n) ~= '/'
  geodir(n+1) = '/';
end;

if abs(L/hbox - round(L/hbox)) > 1e-3 
  disp(['Can not align cube grid with sides, L/hbox = ' num2str(L/hbox) ...
       ' but must be integer']);
  return;
end;
if abs(D/hbox - round(D/hbox)) > 1e-3 
  disp(['Can not align cube grid with sides, D/hbox = ' num2str(D/hbox) ...
       ' but must be integer']);
  return;
end;

disp('Reading geodyn data...');
[t,ur,uz,r,z,rho,cl,ct,pl]=readgeodyndata( geon, geom, geodir, geoname);
disp('Data read. Filter data...');
[tf,vrf,vzf,rf,zf,toff]=filtergeodyn(t,r,z,ur,uz,freq,dt);
disp('Data filtered. Generate whi-file...');
genwhi(filename,tf,rf,zf,vrf,vzf,L,D,hbox,freq,toff);
disp('Whi-file done');
