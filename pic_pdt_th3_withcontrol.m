%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%% 闭环系统状态轨迹演化图 + 2023-10-4 %%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
% lmi_md_th3 % md 判据
lmi_qtd_th2 % qtd 判据
%% 仿真绘图参数
time = 31; % 仿真时长
dt = 0.01; % 仿真步长
N = time/dt; % % 系统迭代次数
sigma = 2; % 快切换区间长度
mode = 2; % 系统模态数
alpha = 0.4; % 指数衰减项
mu = 1.1; % 切换时刻 李雅普诺夫函数的上升率
f = 5;
nu = 1; % 神经元时滞
Ntau = nu/dt;
%% 迭代过程
x_array = zeros(2,N+Ntau); % 系统状态存储矩阵
g_x = zeros(2,N+Ntau); % 激活函数存储矩阵
g_delay_x = zeros(2,N+Ntau); % 延迟 激活函数存储矩阵
w_array = zeros(1,N+Ntau); % 噪声信号存储矩阵
z_array = zeros(1,N+Ntau); % 被控输出存储矩阵
gamma_array = zeros(1,N+Ntau); % 性能指标存储矩阵
%% 系统初值
for i = 1 : Ntau
    x_array(:,i) = [0.3; -0.3];
end
x_array(:,Ntau+1) = [0.3; -0.3];
Mode_array = [];  Phi_array = [];
[Mode,Phi] = ssp_pdt_qtd_mode(time+nu,dt,mode,sigma,alpha,mu,f);
for n = Ntau+1 : N+Ntau
      temp_mode = Mode(n);
      if temp_mode == 1
        A = A1;  B = B1;  C = C1;  W = W1;  tilde_W = tilde_W1;
        H = H1;  E = E1;  D = D1;  L = L1;  J = J1; F = F1;
        eval(strcat('K=K1',int2str(Phi(n)),';'));
%         K = K1;
      else 
        A = A2;  B = B2;  C = C2;  W = W2;  tilde_W = tilde_W2;
        H = H2;  E = E2;  D = D2;  L = L2;  J = J2; F = F2;
        eval(strcat('K=K2',int2str(Phi(n)),';'));
%         K = K2;
      end
      g_x(:,n) = [tanh(x_array(1,n));tanh(x_array(2,n))];
      g_delay_x(:,n) = [tanh(x_array(1,n-Ntau));tanh(x_array(2,n-Ntau))];
      w_array(:,n) = 2.5*exp(-0.25*n*dt)*sin(2.5*pi*n*dt);
      x_array(:,n+1) = x_array(:,n) + dt*( (A + B*K*C + B*H*F*J*C )*x_array(:,n) + W*g_x(:,n) + tilde_W*g_delay_x(:,n) + E*w_array(:,n) );
      z_array(:,n) = D*x_array(:,n) + L*w_array(:,n);
      gamma_array(1,n) = sqrt(norm(z_array(1:n), 2)/norm(w_array(1:n), 2)); % 用于绘制性能指标函数图像
end
disp('**************************************');
disp('************ gamma 实际值 *************');
disp('**************************************'); 
gamma_actual_value = gamma_array(:,N+Ntau)
%% 闭环系统演化轨迹绘制
figure('name', '闭环系统演化轨迹绘制');
plot(0:dt:time-nu-dt, x_array(1, Ntau+1:time/dt), 'color', '#1B5F06', 'linewidth', 1.2, 'linestyle', '-');
hold on
plot(0:dt:time-nu-dt,x_array(2, Ntau+1:time/dt), 'color', '#00AAAA', 'linewidth', 1.2, 'linestyle', '--');
set(gca, 'FontSize', 11, 'Linewidth', 1);
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 15);
ylabel('$x(t)$', 'Interpreter', 'latex', 'fontsize', 15);
h2 = legend('$x_{1}(t)$', '$x_{2}(t)$', 'Interpreter', 'latex', 'fontsize', 15);
set(h2,'Orientation', 'horizon', 'Box', 'on');
grid on;  grid minor;  box on;
axis([0 time-nu, -0.4 0.4]);
%% PDT模式切换信号绘制
figure('name', 'PDT模式切换信号绘制');
plot(0:dt:time-nu-dt, Mode(1, Ntau+1:time/dt), 'color', '#000000', 'linewidth', 1.2);
axis([0, time-nu, 0.5, mode+0.5])
set(gca, 'FontSize', 11, 'Linewidth', 1);
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 15)
ylabel('$\beta(t)$', 'Interpreter', 'latex', 'fontsize', 15)
grid on;  grid minor;  box on;
%% 性能指标绘图
figure('name', '性能指标绘图');
plot(0:dt:time-nu-dt, gamma_array(1, Ntau+1:time/dt), 'color', 'r', 'linewidth', 1.2);   %末尾去掉几个未训练数值
set(gca, 'FontSize', 11, 'Linewidth', 1);
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 15);
ylabel('$\gamma(t)$', 'Interpreter', 'latex', 'fontsize', 15);
legend('$\gamma(t)$', 'Interpreter', 'latex', 'fontsize', 15);
grid on;  grid minor;  box on;
axis([0 time-nu, 0.3 0.48]);
disp(['运行时间: ', num2str(toc)]); 