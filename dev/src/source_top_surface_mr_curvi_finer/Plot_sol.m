% Plot 

clear
clc

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
for jj = 101
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
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(1+5), h1c=contourf(Xci,Xcj,u31cj,clev,'EdgeColor','none');
    hold on
    h1c1=contour(Xci,Xcj,u31cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    h1f = contourf(Xfi,Xfj,u31fj,clev,'EdgeColor','none');
    h1f1 = contour(Xfi,Xfj,u31fj,clev,'LineColor','k'); 
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.1')
    
    clow1 = min(min(u32cj));
    clow2 = min(min(u32fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u32cj));
    chigh2 = max(max(u32fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(2+5), h2c=contourf(Xci,Xcj,u32cj,clev,'EdgeColor','none');
    hold on
    h2c1=contour(Xci,Xcj,u32cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    h2f = contourf(Xfi,Xfj,u32fj,clev,'EdgeColor','none');
    h2f1 = contour(Xfi,Xfj,u32fj,clev,'LineColor','k'); 
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.2')
    
    clow1 = min(min(u33cj));
    clow2 = min(min(u33fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u33cj));
    chigh2 = max(max(u33fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(3+5), h3c=contourf(Xci,Xcj,u33cj,clev,'EdgeColor','none');
    hold on
    h3c1=contour(Xci,Xcj,u33cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    h3f = contourf(Xfi,Xfj,u33fj,clev,'EdgeColor','none');
    h3f1 = contour(Xfi,Xfj,u33fj,clev,'LineColor','k'); 
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.3')
    
    clow1 = min(min(u34cj));
    clow2 = min(min(u34fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u34cj));
    chigh2 = max(max(u34fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(4+5), h4c=contourf(Xci,Xcj,u34cj,clev,'EdgeColor','none');
    hold on
    h4c1=contour(Xci,Xcj,u34cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    h4f = contourf(Xfi,Xfj,u34fj,clev,'EdgeColor','none');
    h4f1 = contour(Xfi,Xfj,u34fj,clev,'LineColor','k'); 
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.4')
    
    clow1 = min(min(u35cj));
    clow2 = min(min(u35fj));
    clow = 0.95*min(clow1,clow2);
    chigh1 = max(max(u35cj));
    chigh2 = max(max(u35fj));
    chigh = 0.95*max(chigh1, chigh2);
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(5+5), h5c=contourf(Xci,Xcj,u35cj,clev,'EdgeColor','none');
    hold on
    h5c1=contour(Xci,Xcj,u35cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    h5f = contourf(Xfi,Xfj,u35fj,clev,'EdgeColor','none');
    h5f1 = contour(Xfi,Xfj,u35fj,clev,'LineColor','k'); 
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.5')
    
    figure(jj)
    h1 = mesh(Xci,Xcj,Xci*0);view(2)
    hold on, h2=mesh(Xfi,Xfj,Xfi*0); view(2)
    set(h1,'EdgeColor','b')
    set(h2,'EdgeColor','r')
    axis equal
    xlabel('x')
    ylabel('z')
    title('y = 1000')
    %print mesh.png
end 


