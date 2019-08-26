% Plot 

clear all

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

u31_c = load('uc3_1.txt'); u31_c = reshape(u31_c,n1_c,n2_c,n3_c);
u31_f = load('uf3_1.txt'); u31_f = reshape(u31_f,n1_f,n2_f,n3_f);
u32_c = load('uc3_2.txt'); u32_c = reshape(u32_c,n1_c,n2_c,n3_c);
u32_f = load('uf3_2.txt'); u32_f = reshape(u32_f,n1_f,n2_f,n3_f);
u33_c = load('uc3_3.txt'); u33_c = reshape(u33_c,n1_c,n2_c,n3_c);
u33_f = load('uf3_3.txt'); u33_f = reshape(u33_f,n1_f,n2_f,n3_f);
u34_c = load('uc3_4.txt'); u34_c = reshape(u34_c,n1_c,n2_c,n3_c);
u34_f = load('uf3_4.txt'); u34_f = reshape(u34_f,n1_f,n2_f,n3_f);
u35_c = load('uc3_5.txt'); u35_c = reshape(u35_c,n1_c,n2_c,n3_c);
u35_f = load('uf3_5.txt'); u35_f = reshape(u35_f,n1_f,n2_f,n3_f);

% this is used to generate 2d mesh (fixed one direction and calculate another two)
for jj = 51
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    u31cj = zeros(n1_c,n3_c);
    u32cj = zeros(n1_c,n3_c);
    u33cj = zeros(n1_c,n3_c);
    u34cj = zeros(n1_c,n3_c); 
    u35cj = zeros(n1_c,n3_c);  
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i);
            Xcj(i,j) = X3c(i,jj,j);
            u31cj(i,j) = u31_c(i,jj,j);
            u32cj(i,j) = u32_c(i,jj,j);
            u33cj(i,j) = u33_c(i,jj,j);
            u34cj(i,j) = u34_c(i,jj,j);
            u35cj(i,j) = u35_c(i,jj,j);
        end 
    end
    Xfi = zeros(n1_f,n3_f);
    Xfj = zeros(n1_f,n3_f);
    u31fj = zeros(n1_f,n3_f);
    u32fj = zeros(n1_f,n3_f);
    u33fj = zeros(n1_f,n3_f);
    u34fj = zeros(n1_f,n3_f);
    u35fj = zeros(n1_f,n3_f);
    for j = 1:n3_f
	      for i = 1:n1_f
            Xfi(i,j) = X1f(i);
            Xfj(i,j) = X3f(i,2*jj-1,j);
            u31fj(i,j) = u31_f(i,2*jj-1,j);
            u32fj(i,j) = u32_f(i,2*jj-1,j);
            u33fj(i,j) = u33_f(i,2*jj-1,j);
            u34fj(i,j) = u34_f(i,2*jj-1,j);
            u35fj(i,j) = u35_f(i,2*jj-1,j);
	      end 
    end
		% find min and max of u31
    clow1 = min(min(u31cj));
    clow2 = min(min(u31fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u31cj));
    chigh2 = max(max(u31fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh,11); % these are the contour levels. You can change the number of levels (11)
    figure(1), h1c=contourf(Xci,Xcj,u31cj,clev);
    axis([0 2000 0 1000])
    hold on
    h1f = contourf(Xfi,Xfj,u31fj,clev);
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.1')
    print t01.png 
    
    % find min and max of u32
    clow1 = min(min(u32cj));
    clow2 = min(min(u32fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u32cj));
    chigh2 = max(max(u32fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh,11); % these are the contour levels. You can change the number of levels (11)
    figure(2), h2c=contourf(Xci,Xcj,u32cj,clev);
    axis([0 2000 0 1000])
    hold on
    h2f = contourf(Xfi,Xfj,u32fj,clev);
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.2')
    print t02.png
    
    % find min and max of u33
    clow1 = min(min(u33cj));
    clow2 = min(min(u33fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u33cj));
    chigh2 = max(max(u33fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh,11); % these are the contour levels. You can change the number of levels (11)
    figure(3), h3c=contourf(Xci,Xcj,u33cj,clev);
    axis([0 2000 0 1000])
    hold on
    h3f = contourf(Xfi,Xfj,u33fj,clev);
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.3')
    print t03.png
    
    % find min and max of u34
    clow1 = min(min(u34cj));
    clow2 = min(min(u34fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u34cj));
    chigh2 = max(max(u34fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh,11); % these are the contour levels. You can change the number of levels (11)
    figure(4), h4c=contourf(Xci,Xcj,u34cj,clev);
    axis([0 2000 0 1000])
    hold on
    h4f = contourf(Xfi,Xfj,u34fj,clev);
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.4')
    print t04.png
    
    % find min and max of u35
    clow1 = min(min(u35cj));
    clow2 = min(min(u35fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u35cj));
    chigh2 = max(max(u35fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh,11); % these are the contour levels. You can change the number of levels (11)
    figure(5), h5c=contourf(Xci,Xcj,u35cj,clev);
    axis([0 2000 0 1000])
    hold on
    h5f = contourf(Xfi,Xfj,u35fj,clev);
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.5')
    print t05.png
    figure(jj)
    h1 = mesh(Xci,Xcj,Xci*0);view(2)
    hold on, h2=mesh(Xfi,Xfj,Xfi*0); view(2)
    set(h1,'EdgeColor','b')
    set(h2,'EdgeColor','r')
    axis equal
    xlabel('x')
    ylabel('z')
    title('y = 1000')
    print mesh.png
end 


