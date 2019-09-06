% Plot 

close all

N6 = load('N6_g.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);
n1_f = N6(4);
n2_f = N6(5);
n3_f = N6(6);

X1c = load('X1c_g.txt'); 
X2c = load('X2c_g.txt'); 
X3c = load('X3c_g.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);

X1f = load('X1f_g.txt'); 
X2f = load('X2f_g.txt'); 
X3f = load('X3f_g.txt'); X3f = reshape(X3f,n1_f,n2_f,n3_f);

err_f = load('err_f_g.txt'); err_f = reshape(err_f,n1_f,n2_f,n3_f);
err_c = load('err_c_g.txt'); err_c = reshape(err_c,n1_c,n2_c,n3_c);

% this is used to generate 2d mesh (fixed one direction and calculate another two)
if (1 == 0)
for jj = 1 
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    errc = zeros(n1_c,n3_c);
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i);
            Xcj(i,j) = X3c(i,jj,j);
		    errc(i,j) = err_c(i,jj,j);
	    end 
    end
    Xfi = zeros(n1_f,n3_f);
    Xfj = zeros(n1_f,n3_f);
    errf = zeros(n1_f,n3_f);
    for j = 1:n3_f
        for i = 1:n1_f
	        Xfi(i,j) = X1f(i);
	        Xfj(i,j) = X3f(i,2*jj-1,j);
		    errf(i,j) = err_f(i,2*jj-1,j);
	    end 
    end
    figure(jj), h1=mesh(Xci,Xcj,Xci*0),view(2)
    hold on, h2=mesh(Xfi,Xfj,Xfi*0), view(2)
	%set(h1,'EdgeColor','b','FaceColor','b','MarkerEdgecolor','b','MarkerFacecolor','b')
	set(h1,'EdgeColor','b')
	set(h2,'EdgeColor','r')
    axis equal
    xlabel('x')
	ylabel('z')
	title('y = 0')
	print -depsc mesh_y0.eps
    %figure(jj+n1_c)
    %mesh(Xci,Xcj,errc)
    %hold on
    %mesh(Xfi,Xfj,errf)
end 
end %if
% this is used to genrate the 3d mesh
if (1 == 0)
[x1c,x2c] = meshgrid(X1c,X2c);
[x1f,x2f] = meshgrid(X1f,X2f);
figure(1)
for i = 1:n3_c
    h3 = mesh(x1c,x2c,X3c(:,:,i));
	set(h3,'EdgeColor','b')
	hold on
end
for i = 1:n3_f
    h4 = mesh(x1f,x2f,X3f(:,:,i));
	set(h4,'EdgeColor','r')
	hold on
end 
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('physical domain')
print -depsc domain.eps
end % if
xc = [0,2,4,6,8];
yc = [0,2,4,6,8];
xf = [0,1,2,3,4,5,6,7,8];
yf = [8,9,10,11,12,13,14,15,16];
xcg = xc;
ycg = [10,10,10,10,10];
[xxc,yyc] = meshgrid(xc,yc);
[xxcg,yycg] = meshgrid(xc,ycg)
[xxf,yyf] = meshgrid(xf,yf);
h6 = mesh(xxf,yyf,0*xxf),view(2)
set(h6,'EdgeColor','r')
hold on
h5 = mesh(xxc,yyc,0*yyc),view(2)
set(h5,'EdgeColor','b')
plot(xcg,ycg,'b*')
%h6 = mesh(xxf,yyf,0*xxf),view(2)
%set(h6,'EdgeColor','r')
grid off
axis off
print -depsc mr.eps

