function [Rx, att, Nx, pe] = IAR_baseline(A, p, lamda, rho)
% 已知基线长度求解基线矢量和整周模糊度
% Integer Ambiguity Resolution
% Rx = [x; y; z]
% att = [psi, theta, rho]基线矢量姿态
% Nx = [N1; N2; ...; Nn]，第一个始终为0
% pe为量测误差的模长，所求的整周模糊度要使它最小
% A的各行为卫星指向天线的单位矢量
% p为两天线相位差不足整周部分，单位：周
% lamda为波长，单位：m
% rho为基线长度，单位：m

% 与第一行做差，排除两天线路径不等长的影响（双差法）
A = A - ones(size(A,1),1)*A(1,:);
A(1,:) = [];
p = p - p(1);
p(1) = [];

% 确定整周模糊度搜索范围
N_max = 2*ceil(rho/lamda); %整周上界，乘2因为上一步的做差
N_min = -N_max; %整周下界

% 搜索
min_pe = inf; %存储当前最小的量测误差
r11 = A(1,1);
r12 = A(1,2);
r13 = A(1,3);
r21 = A(2,1);
r22 = A(2,2);
r23 = A(2,3);
f = r23*r12 - r13*r22;
d2 = (r23*r11 - r13*r21) /  f;
e2 = (r22*r11 - r12*r21) / -f;
a = 1 + d2^2 + e2^2;
for N1=N_min:N_max
    for N2=N_min:N_max
        % 判断方程是否有解
        % r11*x + r12*y + r13*z = s1 = (phi1 + N1)*lamda
        % r21*x + r22*y + r23*z = s2 = (phi2 + N2)*lamda
        % x^2 + y^2 + z^2 = rho^2
        s1 = (p(1)+N1)*lamda;
        s2 = (p(2)+N2)*lamda;
        d1 = (r23*s1 - r13*s2 ) /  f;
        e1 = (r22*s1 - r12*s2 ) / -f;
        b = -2 * (d1*d2 + e1*e2);
        c = d1^2 + e1^2 - rho^2;
        h = b^2-4*a*c;
        if h<0 %方程无解
            continue
        end
        
        for x=[(-b+sqrt(h))/(2*a), (-b-sqrt(h))/(2*a)]
            % 1.计算基线矢量
            y = d1 - d2*x;
            z = e1 - e2*x;
            R = [x;y;z];
            % 2.计算所有整周模糊度
            N = round(A*R/lamda-p);
            % 3.最小二乘计算基线矢量
            R = (A'*A) \ (A'*(p+N)*lamda);
            % 4.比较量测误差
            pe = norm(A*R/lamda-N-p);
            if pe<min_pe
                min_pe = pe;
                Rx = R;
                Nx = N;
            end
        end
    end
end

att = [0,0,0];
att(3) = norm(Rx); %基线长度
att(1) = atan2d(Rx(2),Rx(1)); %基线航向角
att(2) = -asind(Rx(3)/att(3)); %基线俯仰角
Nx = [0; Nx];
pe = min_pe;

end