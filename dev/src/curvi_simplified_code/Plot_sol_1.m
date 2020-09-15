% Plot 

close all

N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);
n1_f = N6(4);
n2_f = N6(5);
n3_f = N6(6);

X1c = load('X1c.txt'); 
X2c = load('X2c.txt'); 
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);

X1f = load('X1f.txt'); 
X2f = load('X2f.txt'); 
X3f = load('X3f.txt'); X3f = reshape(X3f,n1_f,n2_f,n3_f);

% this is used to genrate the 3d mesh
[x2c,x1c] = meshgrid(X1c,X2c);
[x2f,x1f] = meshgrid(X1f,X2f);
figure(1)
for i = 1:n3_c
    h4 = mesh(x1c,x2c,X3c(:,:,i));
	set(h4,'EdgeColor','b','linewidth',1)
	hold on
end

z = zeros(1,n3_c);
for k = 1:n2_c
    y = ones(1,n3_c)*X2c(k);
    for i = 1:n1_c
        x = ones(1,n3_c)*X1c(i);
        for j = 1:n3_c
            z(j) = X3c(i,k,j);
        end 
        hold on
        plot3(x,y,z,'b','linewidth',1)
    end 
end

for i = 1:n3_f
    h4 = mesh(x1f,x2f,X3f(:,:,i));
	set(h4,'EdgeColor','r','linewidth',1)
	hold on
end

z = zeros(1,n3_f);
for k = 1:n2_f
    y = ones(1,n3_f)*X2f(k);
    for i = 1:n1_f
        x = ones(1,n3_f)*X1f(i);
        for j = 1:n3_f
            z(j) = X3f(i,k,j);
        end 
        hold on
        plot3(x,y,z,'r','linewidth',1)
    end 
end

grid off
set(gca,'fontsize',24)
xlabel('x','fontsize',24)
ylabel('y','fontsize',24)
zlabel('z','fontsize',24)

% this is used to generate 2d mesh (fixed one direction and calculate another two)
for jj = 1 
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i);
            Xcj(i,j) = X3c(i,jj,j);
	    end 
    end
    Xfi = zeros(n1_f,n3_f);
    Xfj = zeros(n1_f,n3_f);
    errf = zeros(n1_f,n3_f);
    for j = 1:n3_f
        for i = 1:n1_f
	        Xfi(i,j) = X1f(i);
	        Xfj(i,j) = X3f(i,2*jj-1,j);
	    end 
    end
    figure(2), h1=mesh(Xci,Xcj,Xci*0,'linewidth',1); view(2);
    set(gca,'fontsize',12)
    hold on, h2=mesh(Xfi,Xfj,Xfi*0,'linewidth',1); view(2);
    plot3(X1c,X3f(1:2:end,1,3),X1c(1)*ones(1,n1_c),'b.','MarkerSize',12);view(2)
    set(gca,'fontsize',12)
	set(h1,'EdgeColor','b')
	set(h2,'EdgeColor','r')
    axis equal
    xlabel('x','fontsize',12)
	ylabel('z','fontsize',12)
	%title('y = 0','fontsize',12)
    axis off
	%print -depsc mesh_y0.eps
end 

figure(3)
xc = (0:2:24)*0.5;
yc = (0:4:24)*0.5;
xf = (0:24)*0.5;
yf = (24:2:48)*0.5;
xcg = xc;
ycg = 28*ones(1,13)*0.5;
[xxc,yyc] = meshgrid(xc,yc);
[xxcg,yycg] = meshgrid(xc,ycg);
[xxf,yyf] = meshgrid(xf,yf);
h6 = mesh(xxf,yyf,0*xxf,'linewidth',1);view(2)
set(h6,'EdgeColor','r')
hold on
h5 = mesh(xxc,yyc,0*yyc,'linewidth',1);view(2)
set(h5,'EdgeColor','b')
plot(xcg,ycg,'b.','MarkerSize',12)
axis equal
grid off
axis off