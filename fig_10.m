clc;
clear all;
close all;

tspan = 0:0.01:500;
x0 = rand(1, 3)*49 + 1;
a0 = 2; a1 = 1; b0 = 0.05; c3 = 0.7; w0 = 1; w1= 2; w2 = 1.5; w3 = 2; d0 = 10; d1 = 10; d2 = 10; d3 = 20;
fearo = @(t,v) [a0*v(1) - b0*v(1)^2 - w0*v(1)*v(2)/(d0+v(1));
        w1*v(1)*v(2)/(d1+v(1)) - w2*v(2)*v(3)/(d2+v(2)) - a1*v(2);
        -c3*v(3) + w3*v(2)*v(3)/(d3+v(2))];
[t, v] = ode45(fearo, [0 500], x0);


yp0 = [1; 1; 1];
feari = @(t, y, yp)[yp(1) + 0.1*y(1)*yp(1) - 2.001*y(1) - 0.151*y(1)^2 + 0.005*y(1)^3 + 0.1*y(1)*y(2) ;
   yp(2) + y(2) - 0.099*y(1)*y(2) + 0.15*y(2)*y(3) + 0.099*y(2)^2 + 0.015*y(1)*y(2)*y(3) - 0.01*y(1)*y(2)^2 + 0.099*y(1)*yp(2) + 0.099*y(2)*yp(2) + 0.010*y(1)*y(2)*yp(2);
   yp(3) + 0.7*y(3) - 0.065*y(2)*y(3) + 0.05*y(2)*yp(3)]; 
[x0,yp0] = decic(feari,0,x0,[1 1 0],yp0,[]);
[t,y] = ode15i(feari, [0 500], x0, yp0);
 
figure(1)
[t, v] = ode45(fearo, [0 500], x0);
plot(t,v(:,1),'Linewidth',0.6 )
hold on 
[t,y] = ode15i(feari, [0 500], x0, yp0);
plot(t,y(:,1),'Linewidth',0.6)
legend('Test Data','SINDy-PI model')
xlabel('Time')
ylabel('v_1') 

figure(2)
[t, v] = ode45(fearo, [0 500], x0);
plot(t,v(:,2),'Linewidth',0.6 )
hold on 
[t,y] = ode15i(feari, [0 500], x0, yp0);
plot(t,y(:,2),'Linewidth',0.6)
legend('Test Data','SINDy-PI model');
xlabel('Time')
ylabel('v_2')

figure(3)
[t, v] = ode45(fearo, [0 500], x0);
plot(t,v(:,3),'Linewidth',0.6 )
hold on 
[t,y] = ode15i(feari, [0 500], x0, yp0);
plot(t,y(:,3),'Linewidth',0.6)
legend('Test Data','SINDy-PI model');
xlabel('Time')
ylabel('v_3')

figure(4)
plot3(v(:,1),v(:,2),v(:,3),'Linewidth',0.7 )
hold on 
plot3(y(:,1),y(:,2),y(:,3),'Linewidth',0.7)
legend('Test Data','SINDy-PI model');
xlabel('v_1')
ylabel('v_2')
zlabel('v_3')


