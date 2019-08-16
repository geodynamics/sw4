% Plot 

close all

N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);

X1c = load('X1c.txt'); 
X2c = load('X2c.txt'); 
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);

u31 = load('u3_1.txt'); u31 = reshape(u31,n1_c,n2_c,n3_c);
u32 = load('u3_2.txt'); u32 = reshape(u32,n1_c,n2_c,n3_c);
u33 = load('u3_3.txt'); u33 = reshape(u33,n1_c,n2_c,n3_c);
u34 = load('u3_4.txt'); u34 = reshape(u34,n1_c,n2_c,n3_c);
u35 = load('u3_5.txt'); u35 = reshape(u35,n1_c,n2_c,n3_c);

% this is used to generate 2d mesh (fixed one direction and calculate another two)
if (1 == 1)
for jj = 101
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    u31j = zeros(n1_c,n3_c);
	u32j = zeros(n1_c,n3_c);
	u33j = zeros(n1_c,n3_c);
	u34j = zeros(n1_c,n3_c); 
	u35j = zeros(n1_c,n3_c);  
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i);
            Xcj(i,j) = X3c(i,jj,j);
		    u31j(i,j) = u31(i,jj,j);
			u32j(i,j) = u32(i,jj,j);
			u33j(i,j) = u33(i,jj,j);
			u34j(i,j) = u34(i,jj,j);
			u35j(i,j) = u35(i,jj,j);
	    end 
    end
    figure(1), h1=contourf(Xci,Xcj,u31j)
	colorbar
    xlabel('x')
	ylabel('z')
	title('u3,t=0.1')
	colorbar
	print t01.png 
    figure(2), h2=contourf(Xci,Xcj,u32j)
    colorbar
    xlabel('x')
    ylabel('z')
    title('u3,t=0.2')
	colorbar
	print t02.png
    figure(3), h3=contourf(Xci,Xcj,u33j)
    colorbar
    xlabel('x')
    ylabel('z')
    title('u3,t=0.3')
	colorbar
	print t03.png
    figure(4), h4=contourf(Xci,Xcj,u34j)
    colorbar
    xlabel('x')
    ylabel('z')
    title('u4,t=0.4')
	colorbar
	print  t04.png
    figure(5), h5=contourf(Xci,Xcj,u35j)
    colorbar
    xlabel('x')
    ylabel('z')
    title('u3,t=0.5')
	colorbar
	print  t05.png
end 
end %if

