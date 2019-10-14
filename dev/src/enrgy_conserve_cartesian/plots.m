% Plot 
N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);

X1c = load('X1c.txt'); X1c = reshape(X1c,n1_c,n2_c,n3_c);
X2c = load('X2c.txt'); X2c = reshape(X2c,n1_c,n2_c,n3_c);
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);


uc31 = load('u3_1.txt'); uc31 = reshape(uc31,n1_c,n2_c,n3_c);
jj = 13;
Xci = zeros(n1_c,n3_c);
Xcj = zeros(n1_c,n3_c);
uc31j = zeros(n1_c,n3_c); 
for j = 1:n3_c
    for i = 1:n1_c
        Xci(i,j) = X1c(i,jj,j);
        Xcj(i,j) = X3c(i,jj,j);
		uc31j(i,j) = uc31(i,jj,j);
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
    title('u1,t=0.5')
    colorbar

times = load('times.txt');
energy = load('energy.txt');

figure(2)
plot(times(1:38),energy(1:38),'-o')
