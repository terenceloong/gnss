lamda = 299792458 / 1575.42e6;

baseline = 1; %大致基线长度

% 天线位置
lat = pos(1);
lon = pos(2);
h = pos(3);

p0 = lla2ecef([lat, lon, h]);
Cen = dcmecef2ned(lat, lon);

n = 2000; %计算多少个点，10ms一个点
offset = 0; %偏移多少个点
X = ones(n,3) * NaN; %存结果，第一列航向角，第二列俯仰角，第三列基线长度

for k=1:n %第k个点
    ko = k + offset; %行号
    p = phaseDiffResult(ko,:)'; %相位差
    visibleSV = find(isnan(p)==0); %找到有相位差的通道，通道序号
    svN = length(visibleSV);
    if svN>=4 && ~isnan(posResult(ko,1))%至少四颗卫星
        p = p(~isnan(p)); %删除NaN
        A = zeros(svN,3);
        for m=1:svN %A的第m行
            s = measureResults_A{visibleSV(m)}(ko,1:3); %卫星坐标
            e = (p0-s) / norm(p0-s); %卫星指向天线的方向矢量
            A(m,:) = (Cen*e')'; %转到地理系下
        end
        [Rx, att, Nx, pe] = IAR_nonbaseline(A, p, lamda, baseline+[-0.1,0.1]);
%         [Rx, att, Nx, pe] = IAR_baseline(A, p, lamda, baseline);
        X(k,:) = att;
    end
end