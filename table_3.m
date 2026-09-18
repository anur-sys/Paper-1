clc;
clear all;
close all;
rng(123)

tspan = 0:0.01:1000;
x0 = rand(1, 3)*49 + 1;
x1 = x0 + [0.001 0.001 0.001];
a0 = 2; a1 = 1; b0 = 0.05; c3 = 0.7; w0 = 1; w1= 2; w2 = 1.5; w3 = 2; d0 = 10; d1 = 10; d2 = 10; d3 = 20;
fearo = @(t,v) [a0*v(1) - b0*v(1)^2 - w0*v(1)*v(2)/(d0+v(1));
        w1*v(1)*v(2)/(d1+v(1)) - w2*v(2)*v(3)/(d2+v(2)) - a1*v(2);
        -c3*v(3) + w3*v(2)*v(3)/(d3+v(2))];
[t, v] = ode45(fearo, [0 1000], x0);

% First figure
figure(1)
[t, v1] = ode45(fearo, [0 500], x0);
plot(t, v1(:,1), 'r', 'LineWidth', 0.6)
hold on
[t, v2] = ode45(fearo, [0 500], x1);
plot(t, v2(:,1), 'b', 'LineWidth', 0.6)
legend('(35.1270, 15.0208, 12.1157)', '(35.1280, 15.0218, 12.1167)', 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Time', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_1', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360)
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold

% Second figure
figure(2)
plot3(v(:,1), v(:,2), v(:,3), 'r', 'LineWidth', 0.7)
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold')
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360)
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fear = @(t,v) [2*v(1) - 0.05*v(1)^2 - v(1)*v(2)/(10+v(1));
        1.971*v(1)*v(2)/(10+v(1)) - 1.5*v(2)*v(3)/(10+v(2)) - 0.971*v(2) - 0.294*v(2)/(10+v(1));
        -0.693*v(3) + 1.993*v(2)*v(3)/(20+v(2)) - 0.134*v(3)/(20+v(2))];
[t, v] = ode45(fear, [0 1000], x0);

figure(3)
[t, v1] = ode45(fear, [0 500], x0);
plot(t, v1(:,1),'Linewidth',0.6 )
hold on
[t, v2] = ode45(fear, [0 500], x1);
plot(t,v2(:,1),'Linewidth',0.6 )
legend('(35.1270, 15.0208, 12.1157)','(35.1280, 15.0218, 12.1167)','FontSize', 13, 'FontWeight', 'bold')
xlabel('Time', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_1', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360) 
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold

figure(4)
plot3(v(:,1),v(:,2),v(:,3),'m','Linewidth',0.7 )
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold')
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360)
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yp0 = [1; 1; 1];
feari = @(t, y, yp)[yp(1) + 0.1*y(1)*yp(1) - 2.001*y(1) - 0.151*y(1)^2 + 0.005*y(1)^3 + 0.1*y(1)*y(2) ;
   yp(2) + y(2) - 0.099*y(1)*y(2) + 0.15*y(2)*y(3) + 0.099*y(2)^2 + 0.015*y(1)*y(2)*y(3) - 0.01*y(1)*y(2)^2 + 0.099*y(1)*yp(2) + 0.099*y(2)*yp(2) + 0.010*y(1)*y(2)*yp(2);
   yp(3) + 0.7*y(3) - 0.065*y(2)*y(3) + 0.05*y(2)*yp(3)]; 
[x0,yp0] = decic(feari,0,x0,[1 1 0],yp0,[]);
[t,y] = ode15i(feari, [0 1000], x0, yp0);

figure(5)
[t,y1] = ode15i(feari, [0 500], x0, yp0);
plot(t,y1(:,1),'m','Linewidth',0.6 )
hold on 
[t,y2] = ode15i(feari, [0 500], x1, yp0);
plot(t,y2(:,1), 'b', 'Linewidth',0.6)
legend('(35.1270, 15.0208, 12.1157)','(35.1280, 15.0218, 12.1167)','FontSize', 13, 'FontWeight', 'bold')
xlabel('Time', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_1', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360)
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold

figure(6)
plot3(y(:,1),y(:,2),y(:,3), 'Linewidth',0.7 )
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold')
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360)
set(gca, 'FontSize', 18, 'FontWeight', 'bold') % Adjusts tick size and makes ticks bold
