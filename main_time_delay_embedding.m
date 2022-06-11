% Time delay embedding for Koopman
% from Lecture 3 by J. Nathan Kutz, University of Washington
% YouTube link: https://www.youtube.com/watch?v=K17xYlg_Y_o 

clear all;
close all;

dt = 0.01;
t = 0:dt:50;
x0 = [0.1 5];
mu = 1.2;

rhs = @(t,x)[x(2); mu*(1-x(2)^2)*x(2)-x(1)];

% determine oscillations of a non-linear oscillator by solving an ode
[t,y] = ode45(rhs, t, x0);

plot(t, y(:,1), t, y(:, 2), 'Linewidth', [2]);
xlabel('Time');
ylabel('y1(t), y2(t)');

x1 = y(:, 1);
x2 = y(:, 2);

% Time delay embedding
% define first Henkel matrix with six time delays embedding
H1 = [x1(1:4000).'
      x2(1:4000).'
      x1(2:4001).'
      x2(2:4001).'
      x1(3:4002).'
      x2(3:4002).'
      x1(4:4003).'
      x2(4:4003).'
      x1(5:4004).'
      x2(5:4004).'
      x1(6:4005).'
      x2(6:4005).'];

% define first Henkel matrix with ten time delays embedding
H2 = [x1(1:4000).'
      x2(1:4000).'
      x1(2:4001).'
      x2(2:4001).'
      x1(3:4002).'
      x2(3:4002).'
      x1(4:4003).'
      x2(4:4003).'
      x1(5:4004).'
      x2(5:4004).'
      x1(6:4005).'
      x2(6:4005).'
      x1(7:4006).'
      x2(7:4006).'
      x1(8:4007).'
      x2(8:4007).'
      x1(9:4008).'
      x2(9:4008).'
      x1(10:4009).'
      x2(10:4009).'];

figure(3);
% determine the dominant modes using svd
% there will be 2 dominant modes because it is a 2-D oscillator
% s is the singular value structure
% colums of u are the eigenvectors of the time delays emmbedding Hankel space which is difficult to interpret 
% v is more interpretable than u and conveys information about the dominant correlation structure in time
[u, s, v] = svd(H1, 'econ');
subplot(2,1,1);
plot(diag(s)/(sum(diag(s))), 'o');
xlabel('Time');
ylabel('s');

figure(4);
subplot(2,1,1);
plot(u(:,1:3));
xlabel('Time');
ylabel('u(:,1:3)');
subplot(2,1,2);
plot(v(:, 1:3));
xlabel('Time');
ylabel('v(:, 1:3)');

figure(10), plot(v(:,1:2)); hold on;
xlabel('Time');
ylabel('v(:, 1:3)');

figure(3);
[u, s, v] = svd(H2, 'econ');
subplot(2,1,2), plot(diag(s)/(sum(diag(s))), 'o');

figure(5);
subplot(2,1,1), plot(u(:,1:3));
subplot(2,1,2), plot(v(:,1:3));


H3 = [];
for j = 1:900
    H3 = [H3; y(j:4000+j, :).'];
end


figure(6);
[u,s,v] = svd(H3, 'econ');
plot(diag(s)/(sum(diag(s))), 'o');

figure(7);
subplot(2,1,1), plot(u(:,1:3));
subplot(2,1,2), plot(v(:,1:3));

figure(10), plot(v(:,1:2)); 
legend('6 time delays', '6 time delays', '900 time delays', '900 time delays');





