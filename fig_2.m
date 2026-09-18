clear all;
clc;
close all;

global a0

%a0 = 2;
a1 = 1;
b0 = 0.05;
c3 = 0.7;
w0 = 1;
w1= 2;
w2 = 1.5;
w3 = 2;
d0 = 10;
d1 = 10;
d2 = 10;
d3 = 20;

a0range = 1.5:0.001:2;        % Range for parameter a0
k = 0; 
tspan = 0:0.01:500;          % Time interval for solving food chain system
vmax = [];                  % A matrix for storing the sorted value of v1

for a0 = a0range 
    fear = @(t,v) [a0*v(1) - b0*v(1)^2 - w0*v(1)*v(2)/(d0+v(1));
        w1*v(1)*v(2)/(d1+v(1)) - w2*v(2)*v(3)/(d2+v(2)) - a1*v(2);
        -c3*v(3) + w3*v(2)*v(3)/(d3+v(2))]; 
    v0 = rand(1, 3)*49 + 1;                   
    k = k + 1; 
    [t, v] = ode45(fear, tspan, v0);        % call ode() to solve system
    count = find(t > 100);              % find all the t_values which is > 100 
    v = v(count, :); 
    j = 1; 
    n = length(v(:, 3));                % find the length of vector v1(v in our problem)
    for i = 2:n-1 
        % check for the min value in 1st column of sol matrix
        if (v(i-1, 3) + eps) < v(i, 3) && v(i, 3) > (v(i+1, 3) + eps)
            vmax(k, j) = v(i, 3);       % Sorting the values of v1 in increasing order
            j = j + 1; 
        end 
    end 
    % generating bifurcation map by plotting j-1 element of kth row each time 
    if j > 1 
        plot(a0, vmax(k, 1:j-1), 'b.'); 
    end 
    hold on; 
    index(k) = j-1; 
end 
xlabel('a_{0}'); 
ylabel('{v_{3}}')

% Find transitions and their corresponding h values
threshold = 0.2;
transitions = find(abs(diff(vmax(k,1:j-1))) > threshold);
transition_h = a0range(transitions);

% Add vertical lines to the plot indicating the transitions
%for i = 1:length(transitions)
    %line([transition_h(i) transition_h(i)], [0 max(xmax(k,1:j-1))], 'Color', 'r')
%end

% Print the h values at which transitions occur
disp("Transitions occur at \mu_1 = " + num2str(transition_h))

%Print the number of points plotted for each value of h
disp('Number of points plotted for each value of \mu_1:');
disp(index)

%Print the maximum transient time and the simulation time used for solving food chain system
disp('Maximum transient time:');
disp(max(tspan));
disp('Simulation time:');
disp(length(tspan)*0.1);

% Add plot labels and title
xlabel('a_{0}','FontWeight', 'bold', 'FontSize', 20);
ylabel('{v_{3}}','FontWeight', 'bold', 'FontSize', 20, 'Rotation', 360);
set(gca, 'FontWeight', 'bold', 'FontSize', 16);