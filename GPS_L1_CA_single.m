% 运行双天线前先运行45s单天线，看星历是否正确，还可以预置星历
% 星历不正确会导致码的发射时间不正确，影响定位

%% 计时开始
clear
clc
tic
% 4颗卫星10s数据耗时约16s

%% 文件路径
file_path = 'F:\数据4_30\data_20190430_164934_ch1.dat';
plot_gnss_file(file_path);
sample_offset = 0*4e6;

%% 全局变量
msToProcess = 45*1000; %处理总时间
sampleFreq = 4e6; %接收机采样频率

% 参考位置********************
p0 = [45.730952, 126.624970, 212]; %2A楼顶

buffBlkNum = 40;                     %采样数据缓存块数量（要保证捕获时存储恰好从头开始）
buffBlkSize = 4000;                  %一个块的采样点数（1ms）
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff = zeros(1,buffSize);            %采样数据缓存
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 1.获取文件时间
tf = sscanf(file_path((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')'; %数据文件开始采样时间（日期时间数组）
[tw, ts] = gps_time(tf); %tw：GPS周数，ts：GPS周内秒数
ta = [ts,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间，[s,ms,us]
ta = time_carry(round(ta,2)); %取整

%% 2.根据历书获取当前可能见到的卫星（*）
% svList = [2;6;12];
svList = gps_constellation(tf, p0);
svN = length(svList);

%% 3.为每颗可能见到的卫星分配跟踪通道
channels = repmat(GPS_L1_CA_channel_struct(), svN,1);
for k=1:svN
    channels(k).PRN = svList(k);
    channels(k).state = 0; %状态未激活
end

%% 预置星历
ephemeris_file = ['./ephemeris/',file_path((end-22):(end-8)),'.mat'];
if exist(ephemeris_file, 'file')
    load(ephemeris_file);
    for k=1:svN
        channels(k).ephemeris = ephemeris(k).ephemeris;
    end
end

%% 4.创建跟踪结果存储空间
trackResults = repmat(trackResult_struct(msToProcess+100), svN,1);
for k=1:svN
    trackResults(k).PRN = svList(k);
end

%% 5.创建测量结果存储空间
m = msToProcess/10;
% 接收机时间
receiverTime = zeros(m,3); %[s,ms,us]
% 定位结果
posResult = ones(m,8) * NaN; %位置、速度、钟差、钟频差
%各卫星位置、伪距、速度、伪距率测量结果
measureResults = cell(svN,1);
for k=1:svN
    measureResults{k} = ones(m,8) * NaN;
end

%% 6.打开文件，创建进度条
fclose('all'); %关闭之前打开的所有文件
fileID = fopen(file_path, 'r');
fseek(fileID, round(sample_offset*4), 'bof'); %不取整可能出现文件指针移不过去
if int32(ftell(fileID))~=int32(sample_offset*4)
    error('Sample offset error!');
end
f = waitbar(0, ['0s/',num2str(msToProcess/1000),'s']);

%% 7.信号处理
for t=1:msToProcess
    % 更新进度条
    if mod(t,1000)==0 %1s步进
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    % 读数据（每10s数据用1.2s）
    rawData = double(fread(fileID, [2,buffBlkSize], 'int16')); %取数据，行向量
    buff(buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = rawData(1,:) + rawData(2,:)*1i; %转换成复信号,往缓存里放
    buffBlkPoint = buffBlkPoint + 1;
    buffHead = buffBlkPoint * buffBlkSize;
    if buffBlkPoint==buffBlkNum
        buffBlkPoint = 0; %缓存从头开始
    end
    
    % 更新接收机时间（当前最后一个采样的接收机时间）
% 	ta = time_carry(ta + sample2dt(buffBlkSize, sampleFreq));
    ta = time_carry(ta + [0,1,0]);
    
    %% 捕获（1s搜索一次）
    if mod(t,1000)==0
        for k=1:svN %搜索所有可能见到的卫星
            if channels(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff((end-2*8000+1):end));
                if ~isempty(acqResult) %成功捕获
                    channels(k) = GPS_L1_CA_channel_init(channels(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    trackResults(k).log = [trackResults(k).log; ...
                                             string(['Acquired at ',num2str(t/1000),'s, peakRatio=',num2str(peakRatio)])];
                end
            end
        end
    end
    
    %% 跟踪
    for k=1:svN %第k个通道
        if channels(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults(k).n;
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
                    [channels(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels(k), sampleFreq, buffSize, buff(trackDataTail:trackDataHead));
                else
                    [channels(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels(k), sampleFreq, buffSize, [buff(trackDataTail:end),buff(1:trackDataHead)]);
                end
                % 存跟踪结果（跟踪结果）
                trackResults(k).I_Q(n,:)          = I_Q;
                trackResults(k).disc(n,:)         = disc;
                trackResults(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults(k).CN0(n,:)          = channels(k).CN0;
                trackResults(k).others(n,:)       = others;
                trackResults(k).log               = [trackResults(k).log; log];
                trackResults(k).n                 = n + 1;
            end
        end
    end
    
    %% 定位（每10ms一次）
    if mod(t,10)==0
        n = t/10; %行号
        receiverTime(n,:) = ta; %存接收机时间
        sv = zeros(svN,8);
        for k=1:svN %计算各通道伪距
            if channels(k).state==2 %已经解析到星历
                dn = mod(buffHead-channels(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                codePhase = channels(k).remCodePhase + (dn/sampleFreq)*channels(k).codeFreq; %当前码相位
                ts0 = channels(k).ts0; %码开始的时间，ms
                ts0_code = [floor(ts0/1e3), mod(ts0,1e3), 0]; %码开始的时间，[s,ms,us]，时间用[s,ms,us]表示，避免乘光速时的精度损失
                ts0_phase = [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %码相位时间，[s,ms,us]
                [sv(k,:),~] = sv_ecef(channels(k).ephemeris, ta, ts0_code+ts0_phase); %根据星历计算卫星位置、速度、伪距
                sv(k,8) = -channels(k).carrFreq/1575.42e6*299792458; %载波频率转化为速度
                measureResults{k}(n,:) = sv(k,:); %存储
            end
        end
        sv(sv(:,1)==0,:) = []; %删除没跟踪到的行
        %--------定位--------%
        if size(sv,1)>=4
            pos = pos_solve(sv);
            if abs(pos(7))>0.1 %钟差大于0.1ms时校正接收机钟
                ta = ta - sec2smu(pos(7)/1000);
            else
                posResult(n,:) = pos; %钟准的时候存，时钟不准算的卫星位置有偏差，导致定位不准
            end
        end
    end
    
end

%% 8.关闭文件，关闭进度条
fclose(fileID);
close(f);

%% 删除跟踪结果中的空白数据
for k=1:svN
    trackResults(k) = trackResult_clean(trackResults(k));
end

%% 打印通道日志（*）
clc
for k=1:svN
    if ~isempty(trackResults(k).log)
        disp(['PRN ',num2str(trackResults(k).PRN)])
        n = size(trackResults(k).log,1);
        for kn=1:n
            disp(trackResults(k).log(kn))
        end
        disp(' ')
    end
end
clearvars k n kn

%% 保存星历
ephemeris = struct('PRN',cell(svN,1), 'ephemeris',cell(svN,1));
for k=1:svN
    ephemeris(k).PRN = channels(k).PRN;
    ephemeris(k).ephemeris = channels(k).ephemeris;
end
save(['./ephemeris/',file_path((end-22):(end-8)),'.mat'], 'ephemeris');

%% 画图（*）
for k=1:svN
    if trackResults(k).n==1 %不画没跟踪的通道
        continue
    end
    
    screenSize = get(0,'ScreenSize'); %获取屏幕尺寸
    
    % 建立坐标轴
    %----三图
%     figure('Position', [380, 440, 1160, 480]);
%     ax1 = axes('Position', [0.05, 0.12, 0.4, 0.8]);
%     hold(ax1,'on');
%     axis(ax1, 'equal');
%     title(['PRN = ',num2str(svList(k))])
%     ax2 = axes('Position', [0.5, 0.58, 0.46, 0.34]);
%     hold(ax2,'on');
%     ax3 = axes('Position', [0.5, 0.12, 0.46, 0.34]);
%     hold(ax3,'on');
%     grid(ax3,'on');
    %----五图
    if screenSize(3)==1920 %根据屏幕尺寸设置画图范围
        figure('Position', [390, 280, 1140, 670]);
    elseif screenSize(3)==1368
        figure('Position', [114, 100, 1140, 670]);
    else
        error('Screen size error!')
    end
    ax1 = axes('Position', [0.08, 0.4, 0.38, 0.53]);
    hold(ax1,'on');
    axis(ax1, 'equal');
    title(['PRN = ',num2str(svList(k))])
    ax2 = axes('Position', [0.53, 0.7 , 0.42, 0.25]);
    hold(ax2,'on');
    ax3 = axes('Position', [0.53, 0.38, 0.42, 0.25]);
    hold(ax3,'on');
    grid(ax3,'on');
    ax4 = axes('Position', [0.53, 0.06, 0.42, 0.25]);
    hold(ax4,'on');
    grid(ax4,'on');
    ax5 = axes('Position', [0.05, 0.06, 0.42, 0.25]);
    hold(ax5,'on');
    grid(ax5,'on');
    
    % 画图
    plot(ax1, trackResults(k).I_Q(1001:end,1),trackResults(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.') %I/Q图
    plot(ax2, trackResults(k).dataIndex/sampleFreq, trackResults(k).I_Q(:,1)) %I_P图
    index = find(trackResults(k).CN0~=0);
    plot(ax3, trackResults(k).dataIndex(index)/sampleFreq, trackResults(k).CN0(index), 'LineWidth',2) %载噪比
    
%     index = find(trackResults(k).bitStartFlag==double('H')); %寻找帧头阶段（粉色）
%     plot(ax2, trackResults(k).dataIndex(index)/sampleFreq, trackResults(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     index = find(trackResults(k).bitStartFlag==double('C')); %校验帧头阶段（蓝色）
%     plot(ax2, trackResults(k).dataIndex(index)/sampleFreq, trackResults(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     index = find(trackResults(k).bitStartFlag==double('E')); %解析星历阶段（红色）
%     plot(ax2, trackResults(k).dataIndex(index)/sampleFreq, trackResults(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','r')

    plot(ax4, trackResults(k).dataIndex/sampleFreq, trackResults(k).carrFreq, 'LineWidth',1.5) %载波频率
    plot(ax5, trackResults(k).dataIndex/sampleFreq, trackResults(k).others(:,1)) %视线方向加速度
    
    % 调整坐标轴
    set(ax2, 'XLim',[0,msToProcess/1000])
    
    set(ax3, 'XLim',[0,msToProcess/1000])
    set(ax3, 'YLim',[30,60])
    
    set(ax4, 'XLim',[0,msToProcess/1000])
    set(ax5, 'XLim',[0,msToProcess/1000])
end

clearvars k screenSize ax1 ax2 ax3 ax4 ax5 index

%% 清除变量（*）
clearvars -except sampleFreq msToProcess ...
                  p0 tf svList svN ...
                  channels trackResults ...
                  receiverTime measureResults posResult

%% 计时结束
toc