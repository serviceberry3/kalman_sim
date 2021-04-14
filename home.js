/*

console.log("Welcome");

//Simulate true signal and measurement
//Define system
N = 1000;

//time step
dt = 0.001;

//time vector
t = []
for (i = 0; i < N; i++) {
    t[i] = dt* (i+1);
}*/


/*
F = [1, dt; 0, 1];
G = [-1/2*dt^2; -dt];
H = [1 0];
Q = [0, 0; 0, 0];
u = 9.80665;
I = eye(2);
disp("Test");

y0=100;
v0=0;

xt=zeros(2,N);
xt(:,1)=[y0; v0];

for k=2:N
    xt(:,k)=F*xt(:,k-1)+G*u;
end

R=4;
v=sqrt(R)*randn(1,N);
z=H*xt+v;

x=zeros(2,N);
x(:,1)=[105;0];
P=[10,0;0,0.01];
for k=2:N
    x(:,k)=F*x(:,k-1)+G*u;
    P=F*P*F'+Q;
    K=P*H'/(H*P*H' +R);
    x(:,k)=x(:,k)+K*(z(k)-H*x(:,k));
    P=(I-K*H)*P;
end

figure(1);
subplot(211);
plot(t,z,'g-',t,x(1,:),'b--','LineWidth',2);
hold on; plot(t,xt(1,:),'r:','LineWidth',1.5)
xlabel('t (s)'); ylabel('x_1 = h (m)'); grid on;
legend('Measured', 'Estimated', 'True');
subplot(212);
plot(t, x(2,:), 'b--', 'LineWidth', 2);
hold on; plot(t, xt(2,:), 'r:', 'LineWidth', 1.5)
xlabel('t (s)'); ylabel('x_2 = v (m/s)'); grid on;
legend('Estimated', 'True');

figure(2);
subplot(211);
plot(t, x(1,:)-xt(1,:),'m','LineWidth', 2)
xlabel('t (s)'); ylabel('\Deltax_1 (m)');grid on;
subplot(212);
plot(t, x(2,:)-xt(2,:), 'm', 'LineWidth', 2)
xlabel('t (s)'); ylabel('\Deltax_2 (m/s)'); grid on;*/