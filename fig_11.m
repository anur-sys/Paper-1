clc;
clear all;
close all;
rng(123)

tspan = [0 1000];
x0 = rand(1, 3)*49 + 1;
x1 = x0 + [0.001 0.001 0.001];
a0 = 2; a1 = 1; b0 = 0.05; c3 = 0.7; w0 = 1; w1= 2; w2 = 1.5; w3 = 2; d0 = 10; d1 = 10; d2 = 10; d3 = 20;
fearo = @(t,v) [a0*v(1) - b0*v(1)^2 - w0*v(1)*v(2)/(d0+v(1));
        w1*v(1)*v(2)/(d1+v(1)) - w2*v(2)*v(3)/(d2+v(2)) - a1*v(2);
        -c3*v(3) + w3*v(2)*v(3)/(d3+v(2))];
[t1, x] = ode45(fearo, tspan, x0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fear = @(t,v) [2*v(1) - 0.05*v(1)^2 - v(1)*v(2)/(10+v(1));
        1.971*v(1)*v(2)/(10+v(1)) - 1.5*v(2)*v(3)/(10+v(2)) - 0.971*v(2) - 0.294*v(2)/(10+v(1));
        -0.693*v(3) + 1.993*v(2)*v(3)/(20+v(2)) - 0.134*v(3)/(20+v(2))];
[t2, y] = ode45(fear, tspan, x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yp0 = [1; 1; 1];
feari = @(t, y, yp)[yp(1) + 0.1*y(1)*yp(1) - 2.001*y(1) - 0.151*y(1)^2 + 0.005*y(1)^3 + 0.1*y(1)*y(2) ;
   yp(2) + y(2) - 0.099*y(1)*y(2) + 0.15*y(2)*y(3) + 0.099*y(2)^2 + 0.015*y(1)*y(2)*y(3) - 0.01*y(1)*y(2)^2 + 0.099*y(1)*yp(2) + 0.099*y(2)*yp(2) + 0.010*y(1)*y(2)*yp(2);
   yp(3) + 0.7*y(3) - 0.065*y(2)*y(3) + 0.05*y(2)*yp(3)]; 
[x0,yp0] = decic(feari,0,x0,[1 1 0],yp0,[]);
[t3,z] = ode15i(feari, tspan, x0, yp0);

%Hausdroff metric
function H = hausdorff_distance(X, Y)
    Dxy = pdist2(X, Y);           % distances from X to Y
    Dyx = pdist2(Y, X);           % distances from Y to X

    d_XY = max(min(Dxy, [], 2));  % for each x in X, min dist to any y in Y
    d_YX = max(min(Dyx, [], 2));  % for each y in Y, min dist to any x in X

    H = max(d_XY, d_YX);          % symmetric Hausdorff distance
end

% Example call:
H = hausdorff_distance(x, y);
fprintf('Hausdorff distance between x and y trajectories: %.6f\n', H);

figure(1)
plot3(x(:,1), x(:,2), x(:,3), 'r', 'LineWidth', 0.7);
hold on;
plot3(y(:,1), y(:,2), y(:,3), 'g', 'LineWidth', 0.7);
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360);
legend('Traditional model', 'SINDy model', 'FontSize', 18);


H = hausdorff_distance(x, z);
fprintf('Hausdorff distance between x and z trajectories: %.6f\n', H);

figure(2)
plot3(x(:,1), x(:,2), x(:,3), 'r', 'LineWidth', 0.7);
hold on;
plot3(z(:,1), z(:,2), z(:,3), 'b', 'LineWidth', 0.7);
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360);
legend('Traditional model', 'SINDy-PI model', 'FontSize', 18);

H = hausdorff_distance(y, z);
fprintf('Hausdorff distance between y and z trajectories: %.6f\n', H);
figure(3)
plot3(y(:,1), y(:,2), y(:,3), 'r', 'LineWidth', 0.7);
hold on;
plot3(z(:,1), z(:,2), z(:,3), 'b', 'LineWidth', 0.7);
xlabel('v_1', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('v_2', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('v_3', 'FontSize', 24, 'FontWeight', 'bold','Rotation', 360);
legend('SINDy model', 'SINDy-PI model', 'FontSize', 18);