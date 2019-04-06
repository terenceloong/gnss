%% 计时开始
tic

%% 捕获结果
acqPRN = find(~isnan(acqResults(:,1)))'; %捕获到的卫星编号列表
chN = length(acqPRN); %通道数量

%% 全局变量
msToProcess = 60*1000; %处理总时间
sampleFreq = 4e6; %接收机采样频率

buffBlkNum = 40;                     %采样数据缓存块数量
buffBlkSize = 4000;                  %一个块的采样点数
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff = zeros(1,buffSize);            %采样数据缓存
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 时间
tf = sscanf(file_path((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')';
[~, tf] = gps_time(tf); %数据文件开始采样时间（GPS时间，周内的秒数）
ta = [tf,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间
ta = time_carry(round(ta,2)); %取整

%% 初始化通道
channel = GPS_L1_CA_channel_struct();
channels = repmat(channel, chN,1);
clearvars channel
for k=1:chN
    PRN = acqPRN(k);
    code = GPS_L1_CA_generate(PRN);
    channels(k).PRN = PRN;
    channels(k).n = 1;
    channels(k).state = 0;
    channels(k).trackStage = 'freqPull';
    channels(k).msgStage = 'idle';
    channels(k).cnt = 0;
    channels(k).code = [code(end),code,code(1)];
    channels(k).timeInt = 0.001;
    channels(k).timeIntMs = 1;
    channels(k).codeInt = 1023;
    channels(k).trackDataTail = sampleFreq*0.001 - acqResults(PRN,1) + 2;
    channels(k).blkSize = sampleFreq * 0.001;
    channels(k).trackDataHead = channels(k).trackDataTail + channels(k).blkSize - 1;
    channels(k).dataIndex = channels(k).trackDataTail;
    channels(k).ts0 = NaN;
    channels(k).carrNco = acqResults(PRN,2);
    channels(k).codeNco = 1.023e6 + channels(k).carrNco/1540;
    channels(k).carrAcc = 0;
    channels(k).carrFreq = channels(k).carrNco;
    channels(k).codeFreq = channels(k).codeNco;
    channels(k).remCarrPhase = 0;
    channels(k).remCodePhase = 0;
    channels(k).I_P0 = NaN;
    channels(k).Q_P0 = NaN;
    channels(k).FLL.K = 40;
    channels(k).FLL.Int = channels(k).carrNco;
    [K1, K2] = orderTwoLoopCoef(25, 0.707, 1);
    channels(k).PLL.K1 = K1;
    channels(k).PLL.K2 = K2;
    channels(k).PLL.Int = 0;
    [K1, K2] = orderTwoLoopCoef(2, 0.707, 1);
    channels(k).DLL.K1 = K1;
    channels(k).DLL.K2 = K2;
    channels(k).DLL.Int = channels(k).codeNco;
    channels(k).bitSyncTable = zeros(1,20);
    channels(k).bitBuff = zeros(1,20);
    channels(k).frameBuff = zeros(1,1502);
    channels(k).frameBuffPoint = 0;
    channels(k).ephemeris = zeros(26,1);
    % 计算码鉴相器误差标准差结构体
    channels(k).codeStd.buff = zeros(1,200);
    channels(k).codeStd.buffSize = length(channels(k).codeStd.buff);
    channels(k).codeStd.buffPoint = 0;
    channels(k).codeStd.E0 = 0;
    channels(k).codeStd.D0 = 0;
    % 计算载波鉴相器误差标准差结构体
    channels(k).carrStd.buff = zeros(1,200);
    channels(k).carrStd.buffSize = length(channels(k).carrStd.buff);
    channels(k).carrStd.buffPoint = 0;
    channels(k).carrStd.E0 = 0;
    channels(k).carrStd.D0 = 0;
    channels(k).Px = diag([0.02, 0.01, 5, 1].^2); %6m, 3.6deg, 5Hz, 1Hz/s
end

%% 创建跟踪结果存储空间
dn = 100;
trackResult.PRN = 0;
trackResult.dataIndex     = zeros(msToProcess+dn,1); %码周期开始采样点在原始数据文件中的位置
trackResult.ts0           = zeros(msToProcess+dn,1); %码周期理论发射时间，ms
trackResult.remCodePhase  = zeros(msToProcess+dn,1); %码周期开始采样点的码相位，码片
trackResult.codeFreq      = zeros(msToProcess+dn,1); %码频率
trackResult.remCarrPhase  = zeros(msToProcess+dn,1); %码周期开始采样点的载波相位，周
trackResult.carrFreq      = zeros(msToProcess+dn,1); %载波频率
trackResult.I_Q           = zeros(msToProcess+dn,6); %[I_P,I_E,I_L,Q_P,Q_E,Q_L]
trackResult.disc          = zeros(msToProcess+dn,6); %[codeError,carrError,freqError]
trackResult.bitStartFlag  = zeros(msToProcess+dn,1);
trackResults = repmat(trackResult, chN,1);
clearvars trackResult
for k=1:chN
    trackResults(k).PRN = acqPRN(k);
end

%% 创建测量结果存储空间
measureResults = cell(1,chN+1); %第一列是时间，其他是各个通道
measureResults{1} = zeros(msToProcess/10,3); %接收机时间，[s,ms,us]
for k=1:chN
    measureResults{k+1} = ones(msToProcess/10,8)*NaN;
end
posResult = ones(msToProcess/10,8)*NaN; %存定位结果

%% 打开文件，创建进度条
fileID = fopen(file_path, 'r');
if fileID~=3 %关闭以前打开的文件
    for k=3:(fileID-1)
        fclose(k);
    end
end
fseek(fileID, round(sample_offset*4), 'bof'); %不取整可能出现文件指针移不过去
if int32(ftell(fileID))~=int32(sample_offset*4)
    error('Sample offset error!');
end
f = waitbar(0, ['0s/',num2str(msToProcess/1000),'s']);

%% 信号处理
for t=1:msToProcess
    % 更新进度条
    if mod(t,1000)==0
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    % 1.读数据（每10s数据用1.2s）
    rawData = double(fread(fileID, [2,buffBlkSize], 'int16')); %取数据，行向量
    buff(buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = rawData(1,:) + rawData(2,:)*1i; %转换成复信号,往缓存里放
    buffBlkPoint = buffBlkPoint + 1;
    buffHead = buffBlkPoint * buffBlkSize;
    if buffBlkPoint==buffBlkNum
        buffBlkPoint = 0; %缓存从头开始
    end
    
    % 2.更新接收机时间（当前最后一个采样的接收机时间）
	ta = time_carry(ta + sample2dt(buffBlkSize, sampleFreq));
    
    % 3.通道处理
    for k=1:chN %第k个通道
        while 1
            % 判断是否有完整的跟踪数据
            if mod(buffHead-channels(k).trackDataHead,buffSize)>(buffSize/2)
                break
            end
            
            % 存跟踪结果（通道参数）
            n = channels(k).n;
            trackResults(k).dataIndex(n,:)    = channels(k).dataIndex;
            trackResults(k).ts0(n,:)          = channels(k).ts0;
            trackResults(k).remCodePhase(n,:) = channels(k).remCodePhase;
            trackResults(k).codeFreq(n,:)     = channels(k).codeFreq;
            trackResults(k).remCarrPhase(n,:) = channels(k).remCarrPhase;
            trackResults(k).carrFreq(n,:)     = channels(k).carrFreq;
            
            % 基带处理
            trackDataHead = channels(k).trackDataHead;
            trackDataTail = channels(k).trackDataTail;
            if trackDataHead>trackDataTail
                [channels(k), I_Q, disc, bitStartFlag] = ...
                    GPS_L1_CA_track(channels(k), sampleFreq, buffSize, buff(trackDataTail:trackDataHead));
            else
                [channels(k), I_Q, disc, bitStartFlag] = ...
                    GPS_L1_CA_track(channels(k), sampleFreq, buffSize, [buff(trackDataTail:end),buff(1:trackDataHead)]);
            end
            
            % 存跟踪结果（跟踪结果）
            trackResults(k).I_Q(n,:)          = I_Q; 
            trackResults(k).disc(n,:)         = disc;
            trackResults(k).bitStartFlag(n,:) = bitStartFlag;
        end
    end
    
    % 4.测量结果（每40000个采样点测量一次）
    % 本地钟不准不光影响测距，还会对地球自转坐标变换带来影响，要实时修正钟差和钟频差
    if mod(t,10)==0
        sv = zeros(chN,8);
        n = t/10;
        measureResults{1}(n,:) = ta; %存接收机时间
        for k=1:chN
            if channels(k).state==1 %已经解析到星历
                dn = mod(buffHead-channels(k).trackDataTail+buffSize/2, buffSize) - buffSize/2;
                codePhase = channels(k).remCodePhase + (dn/sampleFreq)*channels(k).codeFreq; %当前码相位
                %----------------------------------------------------------------------------------------------------------%
%                 ts0 = channels(k).ts0/1e3 + codePhase/1.023e6;
%                 tr = ta(1) + ta(2)/1e3 + ta(3)/1e6;
%                 [sv(k,:),~] = sv_ecef_0(channels(k).ephemeris, tr, ts0); %根据星历计算卫星位置、速度、伪距
                %----------------------------------------------------------------------------------------------------------%
                % 时间用[s,ms,us]表示，避免乘光速时的精度损失
                ts0 = channels(k).ts0; %码开始的时间，ms
                ts0_code = [floor(ts0/1e3), mod(ts0,1e3), 0]; %码开始的时间，[s,ms,us]
                ts0_phase = [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %码相位时间，[s,ms,us]
                [sv(k,:),~] = sv_ecef(channels(k).ephemeris, ta, ts0_code+ts0_phase); %根据星历计算卫星位置、速度、伪距
                %----------------------------------------------------------------------------------------------------------%
                sv(k,8) = -channels(k).carrFreq/1575.42e6*299792458; %载波频率转化为速度
                measureResults{k+1}(n,:) = sv(k,:); %存储
            end
        end
        sv(sv(:,1)==0,:) = [];
        if size(sv,1)>=4 %定位
            posResult(n,:) = pos_solve(sv);
        end
    end
end

%% 关闭文件，关闭进度条
fclose(fileID);
close(f);

%% 删除跟踪结果中的空白数据
for k=1:size(trackResults,1)
    n = channels(k).n;
    trackResults(k).dataIndex(n:end,:)    = [];
    trackResults(k).ts0(n:end,:)          = [];
    trackResults(k).remCodePhase(n:end,:) = [];
    trackResults(k).codeFreq(n:end,:)     = [];
    trackResults(k).remCarrPhase(n:end,:) = [];
    trackResults(k).carrFreq(n:end,:)     = [];
    trackResults(k).I_Q(n:end,:)          = [];
    trackResults(k).disc(n:end,:)         = [];
    trackResults(k).bitStartFlag(n:end,:) = [];
end

%% 画I/Q图
for k=1:size(trackResults,1)
    figure
    plot(trackResults(k).I_Q(1001:end,1),trackResults(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.') %画1s之后的图
    axis equal
    title(['PRN = ',num2str(trackResults(k).PRN)])
end
clearvars k

%% 标记比特开始位置
% for k=1:size(trackResults,1)
%     figure
%     plot(trackResults(k).I_Q(:,1))
%     hold on
%     indexBitStart = find(trackResults(k).bitStartFlag==1); %寻找帧头阶段
%     plot(indexBitStart, trackResults(k).I_Q(indexBitStart,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     indexBitStart = find(trackResults(k).bitStartFlag==2); %校验帧头阶段
%     plot(indexBitStart, trackResults(k).I_Q(indexBitStart,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     indexBitStart = find(trackResults(k).bitStartFlag==3); %解析星历阶段
%     plot(indexBitStart, trackResults(k).I_Q(indexBitStart,1), 'LineStyle','none', 'Marker','.', 'Color','r')
%     title(['PRN = ',num2str(trackResults(k).PRN)])
% end
% clearvars k indexBitStart

%% 清除变量
clearvars -except acqResults file_path sample_offset channels trackResults ta measureResults posResult

%% 计时结束
toc