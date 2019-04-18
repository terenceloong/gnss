function test_range()

measureResults = evalin('base', 'measureResults');

% 参考坐标
p0 = lla2ecef([45.74088083, 126.62694533, 197]);
% p0 = lla2ecef([45.741734, 126.62152, 212]);

T = size(measureResults{1},1); %列数
N = length(measureResults); %通道数

d_range = ones(T,N) * NaN; %测距误差
d_vel   = ones(T,N) * NaN; %测速误差
for t=1:T
    for k=1:N
        if ~isnan(measureResults{k}(t,1)) %存在测量结果
            r = measureResults{k}(t,1:3)-p0; %相对位置矢量，接收机指向卫星
            R = norm(r); %真实的距离
            u = r / R; %单位矢量
            dR = dot(measureResults{k}(t,5:7), u); %真实的相对速度，接收机静止

            rho = measureResults{k}(t,4); %测量的距离
            d_range(t,k) = rho - R;

            drho = measureResults{k}(t,8); %测量的速度
            d_vel(t,k) = drho - dR;
        end
    end
end

assignin('base', 'd_range', d_range);
assignin('base', 'd_vel',   d_vel);

figure
plot((1:T)/100, d_range)
legend_text = cell(1,N);
for k=1:N
    legend_text{k} = ['ch', num2str(k)];
end
legend(legend_text)

figure
plot((1:T)/100, d_vel)
legend_text = cell(1,N);
for k=1:N
    legend_text{k} = ['ch', num2str(k)];
end
legend(legend_text)

end