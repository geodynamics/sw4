clear
clc

n1_c = 101;
n2_c = 101;
n3_c = 41;
n1_f = 201;
n2_f = 201;
n3_f = 21;
n1 = 201;
n2 = 201;
n3 = 101;

load('X1c.txt')
load('X2c.txt')
load('X3c.txt')
X3c = reshape(X3c,n1_c,n2_c,n3_c);
load('X1f.txt')
load('X2f.txt')
load('X3f.txt')
X3f = reshape(X3f,n1_f,n2_f,n3_f);

% uniform Cartesian grids (fine)
load('u3_1.txt')
load('u3_2.txt')
load('u3_3.txt')
load('u3_4.txt')
load('u3_5.txt')
u3_1 = reshape(u3_1,n1,n2,n3);
u3_2 = reshape(u3_2,n1,n2,n3);
u3_3 = reshape(u3_3,n1,n2,n3);
u3_4 = reshape(u3_4,n1,n2,n3);
u3_5 = reshape(u3_5,n1,n2,n3);

% curilinear interface with mesh refinement
load('uc3_1.txt')
load('uc3_2.txt')
load('uc3_3.txt')
load('uc3_4.txt')
load('uc3_5.txt')
uc3_1 = reshape(uc3_1,n1_c,n2_c,n3_c);
uc3_2 = reshape(uc3_2,n1_c,n2_c,n3_c);
uc3_3 = reshape(uc3_3,n1_c,n2_c,n3_c);
uc3_4 = reshape(uc3_4,n1_c,n2_c,n3_c);
uc3_5 = reshape(uc3_5,n1_c,n2_c,n3_c);

load('uf3_1.txt')
load('uf3_2.txt')
load('uf3_3.txt')
load('uf3_4.txt')
load('uf3_5.txt')
uf3_1 = reshape(uf3_1,n1_f,n2_f,n3_f);
uf3_2 = reshape(uf3_2,n1_f,n2_f,n3_f);
uf3_3 = reshape(uf3_3,n1_f,n2_f,n3_f);
uf3_4 = reshape(uf3_4,n1_f,n2_f,n3_f);
uf3_5 = reshape(uf3_5,n1_f,n2_f,n3_f);

% calculate difference between two solutions
error1_c = zeros(n1_c,n3_c);
error2_c = zeros(n1_c,n3_c);
error3_c = zeros(n1_c,n3_c);
error4_c = zeros(n1_c,n3_c);
error5_c = zeros(n1_c,n3_c);
% coarse domain
j = 51;
for k = 1:n3_c
    for i = 1:n1_c
        error1_c(i,k) = abs(uc3_1(i,j,k) - u3_1(2*i-1,2*j-1,2*k-1));
        error2_c(i,k) = abs(uc3_2(i,j,k) - u3_2(2*i-1,2*j-1,2*k-1));
        error3_c(i,k) = abs(uc3_3(i,j,k) - u3_3(2*i-1,2*j-1,2*k-1));
        error4_c(i,k) = abs(uc3_4(i,j,k) - u3_4(2*i-1,2*j-1,2*k-1));
        error5_c(i,k) = abs(uc3_5(i,j,k) - u3_5(2*i-1,2*j-1,2*k-1));
    end
end
% fine domain
j = 101;
error1_f = zeros(n1_f,n3_f);
error2_f = zeros(n1_f,n3_f);
error3_f = zeros(n1_f,n3_f);
error4_f = zeros(n1_f,n3_f);
error5_f = zeros(n1_f,n3_f);
for k = 1:n3_f
    for i = 1:n1_f
        error1_f(i,k) = abs(uf3_1(i,j,k) - u3_1(i,j,80+k));
        error2_f(i,k) = abs(uf3_2(i,j,k) - u3_2(i,j,80+k));
        error3_f(i,k) = abs(uf3_3(i,j,k) - u3_3(i,j,80+k));
        error4_f(i,k) = abs(uf3_4(i,j,k) - u3_4(i,j,80+k));
        error5_f(i,k) = abs(uf3_5(i,j,k) - u3_5(i,j,80+k));
    end 
 end
 % l2 error
 l2err1_c = 0; 
 l2err2_c = 0;
 l2err3_c = 0;
 l2err4_c = 0;
 l2err5_c = 0;
for k = 1:n3_c
    for i = 1:n1_c
        l2err1_c = l2err1_c + error1_c(i,k)^2; 
        l2err2_c = l2err2_c + error2_c(i,k)^2; 
        l2err3_c = l2err3_c + error3_c(i,k)^2; 
        l2err4_c = l2err4_c + error4_c(i,k)^2; 
        l2err5_c = l2err5_c + error5_c(i,k)^2; 
     end
 end 
 l2err1_c = sqrt(l2err1_c/n1_c/n3_c);
 l2err2_c = sqrt(l2err2_c/n1_c/n3_c);
 l2err3_c = sqrt(l2err3_c/n1_c/n3_c);
 l2err4_c = sqrt(l2err4_c/n1_c/n3_c);
 l2err5_c = sqrt(l2err5_c/n1_c/n3_c);
 
 l2err1_f = 0; 
 l2err2_f = 0;
 l2err3_f = 0;
 l2err4_f = 0;
 l2err5_f = 0;
for k = 1:n3_f
    for i = 1:n1_f
        l2err1_f = l2err1_f + error1_f(i,k)^2; 
        l2err2_f = l2err2_f + error2_f(i,k)^2; 
        l2err3_f = l2err3_f + error3_f(i,k)^2; 
        l2err4_f = l2err4_f + error4_f(i,k)^2; 
        l2err5_f = l2err5_f + error5_f(i,k)^2; 
     end
 end 
 l2err1_f = sqrt(l2err1_f/n1_f/n3_f);
 l2err2_f = sqrt(l2err2_f/n1_f/n3_f);
 l2err3_f = sqrt(l2err3_f/n1_f/n3_f);
 l2err4_f = sqrt(l2err4_f/n1_f/n3_f);
 l2err5_f = sqrt(l2err5_f/n1_f/n3_f);
 
 % draw the error graphs
 jj = 51;
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
 for j = 1:n3_f
	   for i = 1:n1_f
         Xfi(i,j) = X1f(i);
         Xfj(i,j) = X3f(i,2*jj-1,j);
     end
 end
 figure(1)
 mesh(Xci,Xcj,error1_c), hold on
 mesh(Xfi,Xfj,error1_f)
 xlabel('x')
 ylabel('z')
 zlabel('error')
 title('t01')
 figure(2)
 mesh(Xci,Xcj,error2_c), hold on
 mesh(Xfi,Xfj,error2_f)
 xlabel('x')
 ylabel('z')
 zlabel('error')
 title('t02')
 figure(3)
 mesh(Xci,Xcj,error3_c), hold on
 mesh(Xfi,Xfj,error3_f)
 xlabel('x')
 ylabel('z')
 zlabel('error')
 title('t03')
 figure(4)
 mesh(Xci,Xcj,error4_c), hold on
 mesh(Xfi,Xfj,error4_f)
 xlabel('x')
 ylabel('z')
 zlabel('error')
 title('t04')
 figure(5)
 mesh(Xci,Xcj,error5_c), hold on
 mesh(Xfi,Xfj,error5_f)
 xlabel('x')
 ylabel('z')
 zlabel('error')
 title('t05')
 
 
 