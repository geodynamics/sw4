% Plot 

close all

load('times.txt')
load('energy.txt')

plot(times,energy-energy(1),'-o')