% 递推计算标准差测试程序

X = trackResults(2).disc(:,1);
n = length(X);

buff = zeros(1,200);
buffSzie = length(buff);
buffHead = 0;

result = zeros(n,2);
E0 = 0;
D0 = 0;

for k=1:n
    buffHead = buffHead + 1;
    
    x = X(k);
    x0 = buff(buffHead);
    
    E1 = E0 + (x-x0)/buffSzie;
    D1 = D0 + ((x-E1)^2 - (x0-E0)^2 - 2*(E1-E0)*(E0*buffSzie-x0) + (buffSzie-1)*(E1^2-E0^2))/buffSzie;
    
    result(k,1) = E1;
    result(k,2) = sqrt(D1);
    
    E0 = E1;
    D0 = D1;
    buff(buffHead) = x;
    if buffHead==buffSzie
        buffHead = 0;
    end
end

figure
plot(X)
hold on
plot(result(:,1), 'LineWidth',2)
plot(result(:,2), 'LineWidth',2)