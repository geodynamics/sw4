% Plot 

close all

times = load('times.txt');
energy = load('energy.txt');

plot(times(1:345),energy(1:345),'-o')
