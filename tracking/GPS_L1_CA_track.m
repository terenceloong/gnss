function [channel, I_Q, disc, bitStartFlag] = GPS_L1_CA_track(channel, sampleFreq, buffSize, rawSignal)

%% 载波表
persistent carrTable
if isempty(carrTable)
    carrTable = exp(-(0:3599)/3600*2*pi*1i);
end

%% 输出
bitStartFlag = 0;

%% 更新通道执行次数
channel.n = channel.n + 1;

%% 提取通道信息（跟踪算法用到的控制参数）
trackStage     = channel.trackStage;
msgStage       = channel.msgStage;
cnt            = channel.cnt;
code           = channel.code;
timeInt        = channel.timeInt;
timeIntMs      = channel.timeIntMs;
codeInt        = channel.codeInt;
blkSize        = channel.blkSize;
ts0            = channel.ts0;
carrNco        = channel.carrNco;
codeNco        = channel.codeNco;
carrAcc        = channel.carrAcc;
remCarrPhase   = channel.remCarrPhase;
remCodePhase   = channel.remCodePhase;
I_P0           = channel.I_P0;
Q_P0           = channel.Q_P0;
FLL            = channel.FLL;
PLL            = channel.PLL;
DLL            = channel.DLL;
bitSyncTable   = channel.bitSyncTable;
bitBuff        = channel.bitBuff;
frameBuff      = channel.frameBuff;
frameBuffPoint = channel.frameBuffPoint;
P              = channel.Px; 

%% 基本处理
% 时间序列
t = (0:blkSize) / sampleFreq;
t_2 = t.^2;

% 生成本地载波
% theta = remCarrPhase + carrNco*t;
theta = remCarrPhase + carrNco*t + 0.5*carrAcc*t_2; %变频率载波
% carr = exp(-2*pi*theta(1:end-1)*1i); %本地载波（三角函数）
carr = carrTable(mod(round(theta(1:end-1)*3600),3600)+1); %本地载波（查表）
remCarrPhase = mod(theta(end), 1); %剩余载波相位，周

% 生成本地码
tcode = remCodePhase + codeNco*t;
earlyCode  = code(floor(tcode(1:end-1)+0.5)+2); %超前码
promptCode = code(floor(tcode(1:end-1)    )+2); %即时码
lateCode   = code(floor(tcode(1:end-1)-0.5)+2); %滞后码
remCodePhase = tcode(end) - codeInt; %剩余载波相位，周

% 原始数据乘载波
BasebandSignal = rawSignal .* carr;
iBasebandSignal = real(BasebandSignal);
qBasebandSignal = imag(BasebandSignal);

% 六路积分
I_E = sum(earlyCode  .* iBasebandSignal);
Q_E = sum(earlyCode  .* qBasebandSignal);
I_P = sum(promptCode .* iBasebandSignal);
Q_P = sum(promptCode .* qBasebandSignal);
I_L = sum(lateCode   .* iBasebandSignal);
Q_L = sum(lateCode   .* qBasebandSignal);

% 码鉴相器
S_E = sqrt(I_E^2+Q_E^2);
S_L = sqrt(I_L^2+Q_L^2);
codeError = 0.5 * (S_E-S_L)/(S_E+S_L); %单位：码片
% codeError_coherent = 0.5 * (I_E-I_L)/(I_E+I_L);
[channel.codeStd, codeSigma] = std_rec(channel.codeStd ,codeError); %计算码鉴相器误差标准差

% 载波鉴相器
carrError = atan(Q_P/I_P) / (2*pi); %单位：周
[channel.carrStd, carrSigma] = std_rec(channel.carrStd ,carrError); %计算载波鉴相器误差标准差

% 鉴频器
if ~isnan(I_P0)
    yc = I_P0*I_P + Q_P0*Q_P;
    ys = I_P0*Q_P - Q_P0*I_P;
    freqError = atan(ys/yc)/timeInt / (2*pi); %单位：Hz
else
    freqError = 0;
end

%% 跟踪算法
switch trackStage
    case 'freqPull' %频率牵引
        %----FLL
        FLL.Int = FLL.Int + FLL.K*freqError*timeInt; %锁频环积分器
        carrNco = FLL.Int;
        carrFreq = FLL.Int;
        % 500ms后转到传统跟踪
        cnt = cnt + 1;
        if cnt==500
            cnt = 0; %计数器清零
            PLL.Int = FLL.Int; %初始化锁相环积分器
            trackStage = 'tradTrack';
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
	case 'tradTrack' %传统跟踪
        %----PLL
        PLL.Int = PLL.Int + PLL.K2*carrError*timeInt; %锁相环积分器
        carrNco = PLL.Int + PLL.K1*carrError;
        carrFreq = PLL.Int;
        % 500ms后进行比特同步
        if strcmp(msgStage,'idle')
            cnt = cnt + 1;
            if cnt==500
                cnt = 0; %计数器清零
                msgStage = 'bitSync';
            end
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
    case 'KalmanTrack' %卡尔曼滤波跟踪
        Phi = eye(4);
        Phi(1,3) = timeInt/1540;
        Phi(1,4) = 0.5*timeInt^2/1540;
        Phi(2,3) = timeInt;
        Phi(2,4) = 0.5*timeInt^2;
        Phi(3,4) = timeInt;
        H = zeros(2,4);
        H(1,1) = 1;
        H(2,2) = 1;
        H(2,3) = -timeInt/2;
        Q = diag([0, 0, 3, 0.02].^2) * timeInt^2;
        R = diag([codeSigma, carrSigma].^2);
        Z = [codeError; carrError];
        % 卡尔曼滤波
        P = Phi*P*Phi' + Q;
        K = P*H' / (H*P*H'+R);
        X = K*Z;
        P = (eye(4)-K*H)*P;
        P = (P+P')/2;
        % 修正
        remCodePhase = remCodePhase + X(1);
        remCarrPhase = remCarrPhase + X(2);
        carrNco = carrNco + carrAcc*(blkSize/sampleFreq) + X(3);
        codeNco = 1.023e6 + carrNco/1540;
        carrAcc = carrAcc + X(4);
        carrFreq = carrNco;
        codeFreq = codeNco;
	otherwise
end

%% 电文解析算法
switch msgStage
    case 'bitSync' %比特同步，必须是1ms积分时间，进行2s，有100个比特，比特同步后可以实现更长的积分时间
        cnt = cnt + 1;
        if (I_P0*I_P)<0 %发现电平翻转
            index = mod(cnt-1,20) + 1;
            bitSyncTable(index) = bitSyncTable(index) + 1; %统计表中的对应位加1
        end
        if cnt==2000 %2s后检验统计表
            if max(bitSyncTable)>15 && (sum(bitSyncTable)-max(bitSyncTable))<3 %确定电平翻转位置
                [~,cnt] = max(bitSyncTable);
                cnt = -cnt + 1;
                msgStage = 'findHead';
%                 trackStage = 'KalmanTrack';
            else
                cnt = 0; %计数器清零
                msgStage = 'idle';
            end
            bitSyncTable = zeros(1,20); %比特同步统计表清零
        end
        
    case 'findHead' %寻找帧头
        cnt = cnt + 1; 
        if cnt>0 %往比特缓存中存数
            bitBuff(cnt) = I_P;
        end
        if cnt==1 %标记当前跟踪的数据段为比特开始位置
            bitStartFlag = 1;
        end
        if cnt==(20/timeIntMs) %跟踪完一个比特
            cnt = 0; %计数器清零
            bit = sum(bitBuff(1:(20/timeIntMs))) > 0; %判断比特值，0/1
            frameBuffPoint = frameBuffPoint + 1;
            frameBuff(frameBuffPoint) = (bit - 0.5) * 2; %存储比特值，±1
            %-------------------------------------------------------------%
            if frameBuffPoint>1502 %搜索30s都没有找到帧头
                frameBuffPoint = 0;
                frameBuff = zeros(1,1502); %清空帧缓存
                msgStage = 'idle'; %休息
            elseif frameBuffPoint>=10 %至少有10个比特，前两个用来校验
                if abs(sum(frameBuff(frameBuffPoint+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                    frameBuff(1:10) = frameBuff(frameBuffPoint+(-9:0)); %将帧头提前
                    frameBuffPoint = 10;
                    msgStage = 'checkHead'; %进入校验帧头模式
                end
            end
            %=============================================================%
        end
        
    case 'checkHead' %校验帧头
        cnt = cnt + 1;
        bitBuff(cnt) = I_P; %往比特缓存中存数
        if cnt==1 %标记当前跟踪的数据段为比特开始位置
            bitStartFlag = 2;
        end
        if cnt==(20/timeIntMs) %跟踪完一个比特
            cnt = 0; %计数器清零
            bit = sum(bitBuff(1:(20/timeIntMs))) > 0; %判断比特值，0/1
            frameBuffPoint = frameBuffPoint + 1;
            frameBuff(frameBuffPoint) = (bit - 0.5) * 2; %存储比特值，±1
            %-------------------------------------------------------------%
            if frameBuffPoint==62 %存储了两个字
                if GPS_L1_CA_check(frameBuff(1:32))==1 && GPS_L1_CA_check(frameBuff(31:62))==1 %校验通过
                    % 获取电文时间
                    % bits(32)为上一字的最后一位，校验时控制电平翻转，为1表示翻转，为0表示不翻转，参见ICD-GPS最后几页
                    TOW = -frameBuff(32) * frameBuff(33:49); %31~47比特
                    TOW = bin2dec( dec2bin(TOW>0)' );
                    ts0 = (TOW*6-4.8)*1000 - timeIntMs; %ms
                    msgStage = 'parseEphe'; %进入解析星历模式
                else %校验未通过
                    for ki=11:62 %检查其他比特中有没有帧头
                        if abs(sum(frameBuff(ki+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                            frameBuff(1:10) = frameBuff(ki+(-9:0)); %将帧头提前
                            frameBuffPoint = 10;
                            msgStage = 'checkHead'; %进入校验帧头模式
                        else
                            frameBuff(1:9) = frameBuff(54:62); %将未检测的比特提前
                            frameBuffPoint = 9;
                            msgStage = 'findHead'; %再次寻找帧头
                        end
                    end
                end
            end
            %=============================================================%
        end
        
    case 'parseEphe' %解析星历
        cnt = cnt + 1;
        bitBuff(cnt) = I_P; %往比特缓存中存数
        if cnt==1 %标记当前跟踪的数据段为比特开始位置
            bitStartFlag = 3;
        end
        if cnt==(20/timeIntMs) %跟踪完一个比特
            cnt = 0; %计数器清零
            bit = sum(bitBuff(1:(20/timeIntMs))) > 0; %判断比特值，0/1
            frameBuffPoint = frameBuffPoint + 1;
            frameBuff(frameBuffPoint) = (bit - 0.5) * 2; %存储比特值，±1
            %-------------------------------------------------------------%
            if frameBuffPoint==1502 %跟踪完5帧
                ephemeris = GPS_L1_CA_ephemeris(frameBuff); %解析星历
                if ephemeris(2)==ephemeris(3) %检查星历是否改变
                    channel.ephemeris = ephemeris; %更新星历
                    channel.state = 1; %更新状态
                else
                    disp(['PRN ',num2str(channel.PRN),],' ephemeris change.');
                end
                frameBuff(1:2) = frameBuff(1501:1502); %将最后两个比特提前
                frameBuffPoint = 2;
            end
            %=============================================================%
        end
        
    otherwise
end

%% 更新通道信息1
channel.dataIndex = channel.dataIndex + blkSize;
trackDataHead = channel.trackDataHead;
trackDataTail = trackDataHead + 1;
if trackDataTail>buffSize
    trackDataTail = trackDataTail - buffSize;
end
blkSize = ceil((codeInt-remCodePhase)/codeNco*sampleFreq);
trackDataHead = trackDataTail + blkSize - 1;
if trackDataHead>buffSize
    trackDataHead = trackDataHead - buffSize;
end
channel.ts0           = ts0 + timeIntMs;
channel.trackDataTail = trackDataTail;
channel.blkSize       = blkSize;
channel.trackDataHead = trackDataHead;

%% 更新通道信息2
channel.trackStage     = trackStage;
channel.msgStage       = msgStage;
channel.cnt            = cnt;
channel.code           = code;
channel.timeInt        = timeInt;
channel.timeIntMs      = timeIntMs;
channel.codeInt        = codeInt;
channel.carrNco        = carrNco;
channel.codeNco        = codeNco;
channel.carrAcc        = carrAcc;
channel.carrFreq       = carrFreq;
channel.codeFreq       = codeFreq;
channel.remCarrPhase   = remCarrPhase;
channel.remCodePhase   = remCodePhase;
channel.I_P0           = I_P;
channel.Q_P0           = Q_P;
channel.FLL            = FLL;
channel.PLL            = PLL;
channel.DLL            = DLL;
channel.bitSyncTable   = bitSyncTable;
channel.bitBuff        = bitBuff;
channel.frameBuff      = frameBuff;
channel.frameBuffPoint = frameBuffPoint;
channel.Px             = P;

%% 输出
I_Q = [I_P, I_E, I_L, Q_P, Q_E, Q_L];
% disc = [codeError, carrError, freqError];
% disc = [codeError, codeSigma, carrError, carrSigma, freqError];
disc = [codeError, codeSigma, carrError, carrSigma, freqError, carrAcc];

end