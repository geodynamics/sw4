% Plot 

close all

N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);
n1_f = N6(4);
n2_f = N6(5);
n3_f = N6(6);

X1c = load('X1c.txt'); X1c = reshape(X1c,n1_c,n2_c,n3_c);
X2c = load('X2c.txt'); X2c = reshape(X2c,n1_c,n2_c,n3_c);
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);

X1f = load('X1f.txt'); X1f = reshape(X1f,n1_f,n2_f,n3_f);
X2f = load('X2f.txt'); X2f = reshape(X2f,n1_f,n2_f,n3_f);
X3f = load('X3f.txt'); X3f = reshape(X3f,n1_f,n2_f,n3_f);

err_f = load('err_f.txt'); err_f = reshape(err_f,n1_f,n2_f,n3_f);
err_c = load('err_c.txt'); err_c = reshape(err_c,n1_c,n2_c,n3_c);

for jj = 1:n1_c
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    errc = zeros(n1_c,n3_c);
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i,jj,j);
            Xcj(i,j) = X3c(i,jj,j);
		    errc(i,j) = err_c(i,jj,j);
	    end 
    end
    Xfi = zeros(n1_f,n3_f);
    Xfj = zeros(n1_f,n3_f);
    errf = zeros(n1_f,n3_f);
    for j = 1:n3_f
        for i = 1:n1_f
	        Xfi(i,j) = X1f(i,2*jj-1,j);
	        Xfj(i,j) = X3f(i,2*jj-1,j);
		    errf(i,j) = err_f(i,2*jj-1,j);
	    end 
    end
    figure(jj), mesh(Xci,Xcj,Xci*0),view(2)
    hold on, mesh(Xfi,Xfj,Xfi*0), view(2)
    axis equal
	print mesh.pdf

    figure(jj+n1_c)
    mesh(Xci,Xcj,errc)
    hold on
    mesh(Xfi,Xfj,errf)
end 

