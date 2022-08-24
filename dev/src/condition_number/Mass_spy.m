a = load('Mass.txt');
mass = reshape(a,25*25*3,25*25*3);
figure(1)
spy(mass)
set(gca,'fontsize',12)
figure(2)
spy(mass(1:75,1:75))
set(gca,'fontsize',12)

