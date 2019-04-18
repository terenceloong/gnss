%% 计时开始
clear
clc
tic

%% 文件路径
file_path_A = 'E:\GNSS data\outdoor_static\data_20190410_172036_ch1.dat';
file_path_B = 'E:\GNSS data\outdoor_static\data_20190410_172036_ch2.dat';
sample_offset = 0*4e6;

%% 全局变量
msToProcess = 20*1000; %处理总时间
sampleFreq = 4e6; %接收机采样频率

p0 = [45.74088083, 126.62694533, 197]; %参考位置********************

buffBlkNum = 40;                     %采样数据缓存块数量（要保证捕获时存储恰好从头开始）
buffBlkSize = 4000;                  %一个块的采样点数（1ms）
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff = zeros(2,buffSize);            %采样数据缓存（两行分别为两个天线）
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 1.获取文件时间
tf = sscanf(file_path_A((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')'; %数据文件开始采样时间（日期时间数组）
[tw, ts] = gps_time(tf); %tw：GPS周数，ts：GPS周内秒数
ta = [ts,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间，[s,ms,us]
ta = time_carry(round(ta,2)); %取整

%% 2.根据历书获取当前可能见到的卫星
% svList = [2;6;12];
svList = gps_constellation(tf, p0);
svN = length(svList);

%% 3.为每颗可能见到的卫星分配跟踪通道
channels_A = repmat(GPS_L1_CA_channel_struct(), svN,1);
channels_B = repmat(GPS_L1_CA_channel_struct(), svN,1);
for k=1:svN
    channels_A(k).PRN = svList(k);
    channels_A(k).state = 0; %状态未激活
    channels_B(k).PRN = svList(k);
    channels_B(k).state = 0; %状态未激活
end

%% 预置星历
ephemeris_file = ['./ephemeris/',file_path_A((end-22):(end-8)),'.mat'];
if exist(ephemeris_file, 'file')
    load(ephemeris_file);
    for k=1:svN
        channels_A(k).ephemeris = ephemeris(k).ephemeris;
        channels_B(k).ephemeris = ephemeris(k).ephemeris;
    end
end

%% 4.创建跟踪结果存储空间
trackResults_A = repmat(trackResult_struct(msToProcess+100), svN,1);
trackResults_B = repmat(trackResult_struct(msToProcess+100), svN,1);
for k=1:svN
    trackResults_A(k).PRN = svList(k);
    trackResults_B(k).PRN = svList(k);
end

%% 5.创建测量结果存储空间
m = msToProcess/10;
% 接收机时间
receiverTime = zeros(m,3); %[s,ms,us]
% 定位结果
posResult = ones(m,8) * NaN; %位置、速度、钟差、钟频差
% 相位差
phaseDiffResult = ones(m,svN) * NaN;
%各卫星位置、伪距、速度、伪距率测量结果
measureResults_A = cell(svN,1);
measureResults_B = cell(svN,1);
for k=1:svN
    measureResults_A{k} = ones(m,8) * NaN;
    measureResults_B{k} = ones(m,8) * NaN;
end

%% 6.打开文件，创建进度条
fclose('all'); %关闭之前打开的所有文件
% 文件A
fileID_A = fopen(file_path_A, 'r');
fseek(fileID_A, round(sample_offset*4), 'bof');
if int32(ftell(fileID_A))~=int32(sample_offset*4)
    error('Sample offset error!');
end
% 文件B
fileID_B = fopen(file_path_B, 'r');
fseek(fileID_B, round(sample_offset*4), 'bof');
if int32(ftell(fileID_B))~=int32(sample_offset*4)
    error('Sample offset error!');
end
% 进度条
f = waitbar(0, ['0s/',num2str(msToProcess/1000),'s']);

%% 7.信号处理
for t=1:msToProcess
    % 更新进度条
    if mod(t,1000)==0 %1s步进
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    % 读数据
    rawData = double(fread(fileID_A, [2,buffBlkSize], 'int16'));
    buff(1,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = rawData(1,:) + rawData(2,:)*1i; %天线A
    rawData = double(fread(fileID_B, [2,buffBlkSize], 'int16'));
    buff(2,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = rawData(1,:) + rawData(2,:)*1i; %天线B
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
            %====天线A
            if channels_A(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff(1,(end-2*8000+1):end));
                if ~isempty(acqResult) %成功捕获
                    channels_A(k) = GPS_L1_CA_channel_init(channels_A(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    trackResults_A(k).log = [trackResults_A(k).log; ...
                                             string(['Acquired at ',num2str(t/1000),'s, peakRatio=',num2str(peakRatio)])];
                end
            end
            %====天线B
            if channels_B(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff(2,(end-2*8000+1):end));
                if ~isempty(acqResult) %成功捕获
                    channels_B(k) = GPS_L1_CA_channel_init(channels_B(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    trackResults_B(k).log = [trackResults_B(k).log; ...
                                             string(['Acquired at ',num2str(t/1000),'s, peakRatio=',num2str(peakRatio)])];
                end
            end
        end
    end
    
    %% 跟踪
    for k=1:svN
        %====天线A
        if channels_A(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_A(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_A(k).n;
                trackResults_A(k).dataIndex(n,:)    = channels_A(k).dataIndex;
                trackResults_A(k).ts0(n,:)          = channels_A(k).ts0;
                trackResults_A(k).remCodePhase(n,:) = channels_A(k).remCodePhase;
                trackResults_A(k).codeFreq(n,:)     = channels_A(k).codeFreq;
                trackResults_A(k).remCarrPhase(n,:) = channels_A(k).remCarrPhase;
                trackResults_A(k).carrFreq(n,:)     = channels_A(k).carrFreq;
                % 基带处理
                trackDataHead = channels_A(k).trackDataHead;
                trackDataTail = channels_A(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_A(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq, buffSize, buff(1,trackDataTail:trackDataHead));
                else
                    [channels_A(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq, buffSize, [buff(1,trackDataTail:end),buff(1,1:trackDataHead)]);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_A(k).I_Q(n,:)          = I_Q;
                trackResults_A(k).disc(n,:)         = disc;
                trackResults_A(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_A(k).CN0(n,:)          = channels_A(k).CN0;
                trackResults_A(k).others(n,:)       = others;
                trackResults_A(k).log               = [trackResults_A(k).log; log];
                trackResults_A(k).n                 = n + 1;
            end
        end
        %====天线B
        if channels_B(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_B(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_B(k).n;
                trackResults_B(k).dataIndex(n,:)    = channels_B(k).dataIndex;
                trackResults_B(k).ts0(n,:)          = channels_B(k).ts0;
                trackResults_B(k).remCodePhase(n,:) = channels_B(k).remCodePhase;
                trackResults_B(k).codeFreq(n,:)     = channels_B(k).codeFreq;
                trackResults_B(k).remCarrPhase(n,:) = channels_B(k).remCarrPhase;
                trackResults_B(k).carrFreq(n,:)     = channels_B(k).carrFreq;
                % 基带处理
                trackDataHead = channels_B(k).trackDataHead;
                trackDataTail = channels_B(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_B(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq, buffSize, buff(2,trackDataTail:trackDataHead));
                else
                    [channels_B(k), I_Q, disc, bitStartFlag, others, log] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq, buffSize, [buff(2,trackDataTail:end),buff(2,1:trackDataHead)]);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_B(k).I_Q(n,:)          = I_Q;
                trackResults_B(k).disc(n,:)         = disc;
                trackResults_B(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_B(k).CN0(n,:)          = channels_B(k).CN0;
                trackResults_B(k).others(n,:)       = others;
                trackResults_B(k).log               = [trackResults_B(k).log; log];
                trackResults_B(k).n                 = n + 1;
            end
        end
    end
    
    %% 定位（每10ms一次）暂时使用A天线
    if mod(t,10)==0
        n = t/10; %行号
        receiverTime(n,:) = ta; %存接收机时间
        %--------计算各通道伪距--------%
        sv = zeros(svN,8);
        for k=1:svN
            if channels_A(k).state==2 %已经解析到星历
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                codePhase = channels_A(k).remCodePhase + (dn/sampleFreq)*channels_A(k).codeFreq; %当前码相位
                ts0 = channels_A(k).ts0; %码开始的时间，ms
                ts0_code = [floor(ts0/1e3), mod(ts0,1e3), 0]; %码开始的时间，[s,ms,us]
                ts0_phase = [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %码相位时间，[s,ms,us]
                [sv(k,:),~] = sv_ecef(channels_A(k).ephemeris, ta, ts0_code+ts0_phase); %根据星历计算卫星位置、速度、伪距
                sv(k,8) = -channels_A(k).carrFreq/1575.42e6*299792458; %载波频率转化为速度
                measureResults_A{k}(n,:) = sv(k,:); %存储
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
    
    %% 计算两天线相位差（B-A）（每10ms一次）
    if mod(t,10)==0
        n = t/10; %行号
        for k=1:svN
            if channels_A(k).state==2 && channels_B(k).state==2 %两个天线都跟踪到该颗卫星
                % 天线A
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1; %dn可能等于-1
                dt = dn / sampleFreq;
                phase_A = channels_A(k).remCarrPhase + channels_A(k).carrFreq*dt + 0.5*channels_A(k).carrAcc*dt^2; %载波相位
                % 天线B
                dn = mod(buffHead-channels_B(k).trackDataTail+1, buffSize) - 1; %dn可能等于-1
                dt = dn / sampleFreq;
                phase_B = channels_B(k).remCarrPhase + channels_B(k).carrFreq*dt + 0.5*channels_B(k).carrAcc*dt^2; %载波相位
                % 相位差（小数部分）
                if channels_A(k).inverseFlag*channels_B(k).inverseFlag==1 %两的天线相位翻转相同
                    phaseDiffResult(n,k) = mod(phase_B-phase_A, 1); %单位：周
                else %两的天线相位翻转不同
                    phaseDiffResult(n,k) = mod(phase_B-phase_A+0.5, 1); %单位：周
                end
            end
        end
    end
end

%% 8.关闭文件，关闭进度条
fclose(fileID_A);
fclose(fileID_B);
close(f);

%% 删除跟踪结果中的空白数据
for k=1:svN
    trackResults_A(k) = trackResult_clean(trackResults_A(k));
    trackResults_B(k) = trackResult_clean(trackResults_B(k));
end

%% 打印通道日志
clc
disp('<--------antenna A-------->')
for k=1:svN
    if ~isempty(trackResults_A(k).log)
        disp(['PRN ',num2str(trackResults_A(k).PRN)])
        n = size(trackResults_A(k).log,1);
        for kn=1:n
            disp(trackResults_A(k).log(kn))
        end
        disp(' ')
    end
end
disp('<--------antenna B-------->')
for k=1:svN
    if ~isempty(trackResults_B(k).log)
        disp(['PRN ',num2str(trackResults_B(k).PRN)])
        n = size(trackResults_B(k).log,1);
        for kn=1:n
            disp(trackResults_B(k).log(kn))
        end
        disp(' ')
    end
end

%% 保存星历
ephemeris = struct('PRN',cell(svN,1), 'ephemeris',cell(svN,1));
for k=1:svN
    ephemeris(k).PRN = channels_A(k).PRN;
    if ~isempty(channels_A(k).ephemeris)
        ephemeris(k).ephemeris = channels_A(k).ephemeris;
    else
        ephemeris(k).ephemeris = channels_B(k).ephemeris;
    end
end
save(['./ephemeris/',file_path_A((end-22):(end-8)),'.mat'], 'ephemeris');

%% 画图
for k=1:svN
    if trackResults_A(k).n==1 && trackResults_B(k).n==1 %不画没跟踪的通道
        continue
    end
    
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
    %----五图
    figure('Position', [390, 280, 1140, 670]);
    ax1 = axes('Position', [0.08, 0.4, 0.38, 0.53]);
    hold(ax1,'on');
    axis(ax1, 'equal');
    title(['PRN = ',num2str(svList(k))])
    ax2 = axes('Position', [0.53, 0.7 , 0.42, 0.25]);
    hold(ax2,'on');
    ax3 = axes('Position', [0.53, 0.38, 0.42, 0.25]);
    hold(ax3,'on');
    ax4 = axes('Position', [0.53, 0.06, 0.42, 0.25]);
    hold(ax4,'on');
    grid(ax4,'on');
    ax5 = axes('Position', [0.05, 0.06, 0.42, 0.25]);
    hold(ax5,'on');
    grid(ax5,'on');
    
    % 画图
    plot(ax1, trackResults_A(k).I_Q(1001:end,1),trackResults_A(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.', 'Color',[0,0.447,0.741])
    plot(ax2, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).I_Q(:,1), 'Color',[0,0.447,0.741])
    
%     index = find(trackResults_A(k).bitStartFlag==double('H')); %寻找帧头阶段（粉色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     index = find(trackResults_A(k).bitStartFlag==double('C')); %校验帧头阶段（蓝色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     index = find(trackResults_A(k).bitStartFlag==double('E')); %解析星历阶段（红色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','r')
    %---------------------------------------------------------------------%
    plot(ax1, trackResults_B(k).I_Q(1001:end,1),trackResults_B(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.', 'Color',[0.850,0.325,0.098])
    plot(ax3, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).I_Q(:,1), 'Color',[0.850,0.325,0.098])
    
%     index = find(trackResults_B(k).bitStartFlag==double('H')); %寻找帧头阶段（粉色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     index = find(trackResults_B(k).bitStartFlag==double('C')); %校验帧头阶段（蓝色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     index = find(trackResults_B(k).bitStartFlag==double('E')); %解析星历阶段（红色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','r')

    plot(ax4, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).carrFreq, 'LineWidth',1.5, 'Color',[0,0.447,0.741]) %载波频率
    plot(ax4, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).carrFreq, 'LineWidth',1.5, 'Color',[0.850,0.325,0.098])
    
    plot(ax5, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).others(:,1), 'Color',[0,0.447,0.741]) %视线方向加速度
    plot(ax5, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).others(:,1), 'Color',[0.850,0.325,0.098])
    
    % 调整坐标轴
    set(ax2, 'XLim',[0,msToProcess/1000])
    set(ax3, 'XLim',[0,msToProcess/1000])

    ax2_ylim = get(ax2, 'YLim');
    ax3_ylim = get(ax3, 'YLim');
    ylim = max(abs([ax2_ylim,ax3_ylim]));
    set(ax2, 'YLim',[-ylim,ylim])
    set(ax3, 'YLim',[-ylim,ylim])
    
    set(ax4, 'XLim',[0,msToProcess/1000])
    set(ax5, 'XLim',[0,msToProcess/1000])
end

%% 计时结束
toc