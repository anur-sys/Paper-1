clc;
clear all;
close all;

tspan = 0:0.01:1000;
x0 = [35.23895, 10.76923, 7.72481];
x1 = [35.23905, 10.76933, 7.72491];
a0 = 2; a1 = 1; b0 = 0.05; c3 = 0.7; w0 = 1; w1= 2; w2 = 1.5; w3 = 2; d0 = 10; d1 = 10; d2 = 10; d3 = 20;
fearo = @(t,v) [a0*v(1) - b0*v(1)^2 - w0*v(1)*v(2)/(d0+v(1));
        w1*v(1)*v(2)/(d1+v(1)) - w2*v(2)*v(3)/(d2+v(2)) - a1*v(2);
        -c3*v(3) + w3*v(2)*v(3)/(d3+v(2))];
[t, v] = ode45(fearo, [0 500], x0);

figure(1)
[t, v] = ode45(fearo, [0 700], x0);
plot(t,v(:,1),'Linewidth',0.6 )
hold on
[t, v1] = ode45(fearo, [0 700], x1);
plot(t,v1(:,1),'Linewidth',0.6 )
legend('(35.23895, 10.76923, 7.72481)','(35.23905, 10.76933, 7.72491)')
xlabel('Time')
ylabel('v_1') 

figure(2)
plot3(v(:,1),v(:,2),v(:,3),'Linewidth',0.7 )
xlabel('v_1')
ylabel('v_2')
zlabel('v_3')