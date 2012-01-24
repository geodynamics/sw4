% function loh3exact
%
% usage: [t,ra,tr,ve] = loh3exact(sig)
%
%  Input:
%    sig:  spread in Gaussian time function (0.05 for original LOH3 test)
% 
%  Output:
%    t:    vector of time levels
%    ra:   vector of radial velocities
%    tr:   vector of transverse velocities
%    ve:   vector of vertical velocities (NOTE: positive is down)
%
    function [t,ra,tr,ve] = loh3exact(sig)
filename = 'LOH.3_prose_corrected';
fid=fopen(filename);
A=fscanf(fid,'%g',inf);
nt=length(A)/4;
B=reshape(A,4,nt);
t=B(1,:);
dt=t(2)-t(1);
vert1=B(2,:)*(-1e5);
rad1=B(3,:)*1e5;
trans1=B(4,:)*1e5;
fclose(fid);
%rad=real(ifft(conj(fft(rad(1:nt)))));
%trans=real(ifft(conj(fft(trans(1:nt)))));
%vert=real(ifft(conj(fft(vert(1:nt)))));
T=0.1;
%
filter=(1/T^2)*t.*exp(-t./T);
ra=dt*conv(rad1,filter);
tr=dt*conv(trans1,filter);
ve=dt*conv(vert1,filter);
ve=ve(1:nt);
ra=ra(1:nt);
tr=tr(1:nt);
T=0.1;ts=4*sig;tau=t-ts;
factor=1-(2*T/sig^2)*tau-((T/sig)^2)*(1-(tau./sig).^2);
filter=(1/sqrt(2*pi)/sig)*factor.*exp(-0.5*(tau./sig).^2);
ra=dt*conv(ra,filter);
tr=dt*conv(tr,filter);
ve=dt*conv(ve,filter);
%
% NOTE NO CHANGE OF SIGN!
ve=ve(1:nt); 
ra=ra(1:nt);
tr=tr(1:nt);

