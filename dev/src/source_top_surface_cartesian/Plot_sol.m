% Plot 

clear 
clc

N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);

X1c = load('X1c.txt'); 
X2c = load('X2c.txt'); 
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);


uc31 = load('u2_1.txt'); uc31 = reshape(uc31,n1_c,n2_c,n3_c);
uc32 = load('u2_2.txt'); uc32 = reshape(uc32,n1_c,n2_c,n3_c);
uc33 = load('u2_3.txt'); uc33 = reshape(uc33,n1_c,n2_c,n3_c);
uc34 = load('u2_4.txt'); uc34 = reshape(uc34,n1_c,n2_c,n3_c);
uc35 = load('u2_5.txt'); uc35 = reshape(uc35,n1_c,n2_c,n3_c);


% this is used to generate 2d mesh (fixed one direction and calculate another two)
jj = 101;
Xci = zeros(n1_c,n3_c);
Xcj = zeros(n1_c,n3_c);
uc31j = zeros(n1_c,n3_c);
uc32j = zeros(n1_c,n3_c);
uc33j = zeros(n1_c,n3_c);
uc34j = zeros(n1_c,n3_c); 
uc35j = zeros(n1_c,n3_c);  
for j = 1:n3_c
    for i = 1:n1_c
        Xci(i,j) = X1c(i);
        Xcj(i,j) = X3c(i,jj,j);
		    uc31j(i,j) = uc31(i,jj,j);
			  uc32j(i,j) = uc32(i,jj,j);
			  uc33j(i,j) = uc33(i,jj,j);
			  uc34j(i,j) = uc34(i,jj,j);
			  uc35j(i,j) = uc35(i,jj,j);
	  end 
end

    % find min and max of u31
    clow1 = min(min(uc31j));
    clow = 0.95*clow1;
    chigh1 = max(max(uc31j));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh,10); % these are the contour levels. You can change the number of levels (11)
    figure(1), h1c=contourf(Xci,Xcj,uc31j,clev);
    axis equal
    xlabel('x')
    ylabel('z')
    title('u2,t=0.1')
    colorbar
    
    % find min and max of u32
    clow1 = min(min(uc32j));
    clow = 0.95*clow1;
    chigh1 = max(max(uc32j));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh,10); % these are the contour levels. You can change the number of levels (11)
    figure(2), h2c=contourf(Xci,Xcj,uc32j,clev);
    axis equal
    xlabel('x')
    ylabel('z')
    title('u2,t=0.2')
    colorbar
    
    % find min and max of u33
    clow1 = min(min(uc33j));
    clow = 0.95*clow1;
    chigh1 = max(max(uc33j));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh,10); % these are the contour levels. You can change the number of levels (11)
    figure(3), h3=contourf(Xci,Xcj,uc33j,clev);
    axis equal
    xlabel('x')
    ylabel('z')
    title('u2,t=0.3')
    colorbar
    
    % find min and max of u34
    clow1 = min(min(uc34j));
    clow = 0.95*clow1;
    chigh1 = max(max(uc34j));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh,10); % these are the contour levels. You can change the number of levels (11)
    figure(4), h4=contourf(Xci,Xcj,uc34j,clev);
    axis equal
    xlabel('x')
    ylabel('z')
    title('u2,t=0.4')
    colorbar
    
    % find min and max of u35
    clow1 = min(min(uc35j));
    clow = 0.95*clow1;
    chigh1 = max(max(uc35j));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh,10); % these are the contour levels. You can change the number of levels (11)
    figure(5), h5=contourf(Xci,Xcj,uc35j,clev);
    axis equal
    xlabel('x')
    ylabel('z')
    title('u2,t=0.5')
    colorbar
    
    figure(jj)
    h1 = mesh(Xci,Xcj,Xci*0);view(2)
    set(h1,'EdgeColor','r')
    axis equal
    xlabel('x')
    ylabel('z')
    title('y = 1000')

