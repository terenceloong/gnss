% 测试短基线整周模糊度搜索算法
% 观察R和Rx是否相等、att和设置的姿态是否相等、N和Nx差同一个常数

%% 生成数据
lamda = 299792458 / 1575.42e6;

% 卫星的方位角、高度角
sv = [0,   45;
      23,  28;
      58,  80;
      100, 49;
      146, 34;
      186, 78;
      213, 43;
      255, 86;
      310, 20];

% 基线矢量
psi = 270;
theta = 22;
rho = 1.9;
R = [cosd(theta)*cosd(psi); cosd(theta)*sind(psi); -sind(theta)] * rho;

n = size(sv,1);
A = zeros(n,3);
p = zeros(n,1);
N = zeros(n,1);
for k=1:n
    A(k,:) = [-cosd(sv(k,2))*cosd(sv(k,1)); -cosd(sv(k,2))*sind(sv(k,1)); sind(sv(k,2))]; %卫星指向天线
    phase = dot(A(k,:),R) / lamda + randn(1)*0.005; %相位差加噪声，单位：周
    N(k) = floor(phase); %整数部分
    p(k) = mod(phase,1); %小数部分
end

%% 求解整周模糊度
% [att, Rx, Nx, pe] = IAR_baseline(A, p, lamda, rho); %已知基线求解整周模糊度
[att, Rx, Nx, i1, pe] = IAR_nonbaseline(A, p, lamda, rho+[-0.1,0.1]); %未知基线求解整周模糊度

if length(unique(N-Nx))~=1 %N-Nx应该都相同，如果不同说明算错了
    error('Error!')
end