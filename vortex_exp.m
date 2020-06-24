close all;clear;clc;
format long

% ---------------------------------------- 控制台  ----------------------------------------
        RLS = 2;                                                           % 0 --无均衡   | 1 --RLS均衡  | 2 --LS估计迫零均衡
        item = 1;                                                          % 实验步骤(1--case1; 2~4--case2)
        status = 2;                                                        % 0 --测试状态 | 1 --实验状态
        Fs = 50e6;                                                         % 采样频率
        Cycles = 1;                                                        % 重复次数
        TimeSlot = 30;                                                     % 数据总量
        
        sigma = 0.6;                                                      % RLS 参数[0.8:2]
        Order =3;
        
        Cut = 0;                                                           % 剪切标志位
        factor = 13.5;                                                     % 剪切率
        
        Folder = '用户实验';                                                % 存放目录
        OveSamT = 4;                                                       % 过采样倍数
        OffLineData = '用户实验\ExpData1.mat';                            % 离线数据目录
        
% --------------------------------------------------------------------------------------------

modebit = containers.Map;
curmode = 'QPSK';
modebit('QPSK')  = 2; 
modebit('16QAM') = 4; 
modebit('32QAM') = 5;
modebit('64QAM') = 6; 
modebit('128QAM') = 7; 
modebit('256QAM') = 8;

PrefixRatio = 1/8;
SymbolNum = 24;                                                            % OFDM 符号数
IFFTLength = 128;
GI = PrefixRatio*IFFTLength;
SubCarrier = IFFTLength/2-1;                                               % 最大化利用子载波
BitSymbol  = modebit(curmode);
EvmRate = zeros(1,length(Cycles));                                         % EVM
EVMsub = zeros(Cycles,SubCarrier);
EVMuser= zeros(Cycles,3);
SymbolLen = IFFTLength+GI;
TimeSlotLen = (SymbolNum+1)*SymbolLen;
FrameLon = TimeSlot*SubCarrier*SymbolNum*BitSymbol;                        % 总数据量
UserTotal = TimeSlot*SubCarrier*SymbolNum/3;
disp(['Total simulation: ',num2str(FrameLon)]);
disp(' ');

% 同步序列
N = 63;
n = 1:1:N;
nzc = exp( -1i*pi*34*(n-1).*n./N);                                         % ZC序列
Nzc = [nzc,nzc(1)];
NzcJ = conj(Nzc);
SynSeq = zeros(1,IFFTLength);
SynSeq(2:1+length(Nzc)) = Nzc;
SynSeq(end:-1:end-length(Nzc)+1) = NzcJ(1:end);                            % 共轭部分
SynBlock = ifft(SynSeq);
SynPart = [SynBlock(end-GI+1:end),SynBlock];                               % 添加CP

% 导频
N = SubCarrier-1;
n = 1:1:N;
nzc = exp( -1i*pi*25*(n-1).*n./N);
Nzc = [nzc,nzc(1)];
PilotSeq = zeros(1,IFFTLength);
PilotSeq(2:1+length(Nzc)) = Nzc;
NzcJ = conj(Nzc);
PilotSeq(end:-1:end-length(Nzc)+1) = NzcJ(1:end);   % 共轭
PilotBlock = ifft(PilotSeq);
% -------------------------------------------------------  发送端  ------------------------------------------------------- 

if exist('用户实验','dir') == 0
   mkdir('用户实验'); 
end

if exist('用户实验\Info.txt','file') == 0  
    Num = 0;
    FileName = sprintf('%s\\Info.txt',Folder);
    save(FileName, 'Num', '-ASCII', '-double');
    
    TranData = randi([0,1],1,FrameLon);
    DataName = sprintf('%s\\TranData.txt',Folder);
    save(DataName, 'TranData', '-ASCII', '-double');
    MapDataT = modulate_def(TranData,curmode,1);
    FrameMap = reshape(MapDataT,SubCarrier,TimeSlot*SymbolNum);
    FrameMap = FrameMap.';
    
    for user =1:3
        if user==1
            FrameMapT = FrameMap;
            FrameMapT(:,22:63) = 0;
            User1MapT = FrameMapT(:,1:21);            
            UserMap = reshape(User1MapT.',1,FrameLon/3/BitSymbol);      
        elseif user==2
            FrameMapT = FrameMap;
            FrameMapT(:,1:21) = 0;
            FrameMapT(:,43:63) = 0;
            User2MapT = FrameMapT(:,22:42);
            UserMap = reshape(User2MapT.',1,FrameLon/3/BitSymbol);
        elseif user==3
            FrameMapT = FrameMap;
            FrameMapT(:,1:42) = 0;
            User3MapT = FrameMapT(:,43:63);            
            UserMap = reshape(User3MapT.',1,FrameLon/3/BitSymbol);
        end
        
        % 时域变换
        FrameBox = zeros(TimeSlot*SymbolNum,IFFTLength);
        FrameBox(:,2:SubCarrier+1) = FrameMapT;
        FrameBox(:,end:-1:end-SubCarrier+1) = conj(FrameMapT);
        FrameBox = ifft(FrameBox,IFFTLength,2);
        FrameTran = zeros(TimeSlot*(SymbolNum+1),IFFTLength);
            
        Peak = max(max(abs(FrameBox)));                                    % 归一化峰值，消除突出毛刺
        SynPart = Peak*SynPart./max(SynPart)./3;
        PilotBlock = Peak*PilotBlock./max(PilotBlock)./3;
                
        % 导频及前缀
        for slot = 0:TimeSlot-1
            FrameTran(slot*(SymbolNum+1)+1,:) = PilotBlock;
            FrameTran((slot*(SymbolNum+1)+2):(slot*(SymbolNum+1)+1+SymbolNum),:) = FrameBox(slot*SymbolNum+1:(slot+1)*SymbolNum,:);
        end
        FramePort = [FrameTran(:,IFFTLength-GI+1:end),FrameTran];
        FramePort = reshape(FramePort.',1,TimeSlotLen*TimeSlot);

        % 同步序列
        TranLine = [SynPart,FramePort];               
        SignalBackUp = TranLine;

        % 过采样[4 times]
        Len = length(TranLine);
        TranLine = oversamp(TranLine, Fs, Fs*OveSamT, 0, 0,(OveSamT - 1)*Len);
        TranLine = real(TranLine);

        % 信号剪切  
        if Cut == 1
            factor = 10^(factor/10);
            signalpower = mean(TranLine.^2);
            ClippingThreshold = sqrt(factor*signalpower);                  % 求剪切门限
            TranLine(TranLine> ClippingThreshold) = ClippingThreshold;
            TranLine(TranLine<-ClippingThreshold) =-ClippingThreshold;
        end
                
        % 存储发送信号
        TranLine = TranLine.';
        DataName = sprintf('%s\\User%d Signal.txt',Folder,user);
        save(DataName, 'TranLine', '-ASCII', '-double');
    end
    
    % 存储星座图
    FrameMapr = real(FrameMap);
    FrameMapi = imag(FrameMap);
    MapName = sprintf('%s\\FrameMapr.txt',Folder);
    save(MapName, 'FrameMapr', '-ASCII', '-double');
    MapName = sprintf('%s\\FrameMapi.txt',Folder);
    save(MapName, 'FrameMapi', '-ASCII', '-double');
    DataLen = length(TranLine);
    disp(['生成数据量为：',num2str(DataLen)]);
end


% ---------------------------------------- 加载离线数据  ----------------------------------------
DataName = sprintf('%s\\FrameMapr.txt',Folder);
FrameMapr = load(DataName);
DataName = sprintf('%s\\FrameMapi.txt',Folder);
FrameMapi = load(DataName);
FrameMapT  = FrameMapr + 1i.*FrameMapi;

User1MapT = FrameMapT(:,1:21);            
User1MapT = reshape(User1MapT.',1,FrameLon/3/BitSymbol);

User2MapT = FrameMapT(:,22:42);            
User2MapT = reshape(User2MapT.',1,FrameLon/3/BitSymbol);

User3MapT = FrameMapT(:,43:63);            
User3MapT = reshape(User3MapT.',1,FrameLon/3/BitSymbol);
    
if status == 0
    disp('----------------------当前为测试状态----------------------')
    DataName = sprintf('%s\\User1 Signal.txt',Folder);
    UserSignal1 = importdata(DataName);
    
    DataName = sprintf('%s\\User2 Signal.txt',Folder);
    UserSignal2 = importdata(DataName);
    
    DataName = sprintf('%s\\User3 Signal.txt',Folder);
    UserSignal3 = importdata(DataName);
    
    disp('数据导入完毕!');
    TranLine = UserSignal1 + UserSignal2 + UserSignal3;
    TranLine = TranLine.';
elseif status == 1                                                         % 数据导入
    disp('----------------------当前为实验状态----------------------')
    Point = 2e6;
    OSC = OSCInit();        % 连接示波器
    OSCPoint(OSC,Point);    % 设置采样点
elseif status == 2
    disp('----------------------当前为离线状态----------------------')
end

for cyc = 1:Cycles
    disp(['当前循环:',num2str(cyc)]);
    if status == 1                                                         % 实验状态
        TranLine = OSCRead(OSC);
        disp(['长度：',num2str(length(TranLine))]);
        
        DataName = sprintf('%s\\Info.txt',Folder);
        Num = load(DataName);
        Num = Num + 1;
        filename = sprintf('%s\\ExpData%d',Folder,Num);
        save(strcat(filename,'.mat'),'TranLine');
        save(DataName, 'Num', '-ASCII', '-double');
        disp('收集数据存储完毕！');
    elseif status == 2                                                     % 离线状态
%         OffLineData = sprintf('用户实验\\ExpData%d.mat',cyc);
        as = load(OffLineData);
        TranLine = as.TranLine;
    end

    % ----------------------------------------------      LPF     ----------------------------------------------
    lpf = zeros(1,length(TranLine));
    lpf(1:1:ceil(end/2/OveSamT)) = 1;
    lpf(end:-1:ceil(end-end/2/OveSamT+1)) = 1;
    TranLine = real(ifft(lpf.*fft(TranLine)));
    
%     TranLine = LMS(TranLine);
    % 过采样
    RecLine = TranLine(1:OveSamT:end);
    
    RecLine = RecLine - mean(RecLine);
    [FrameLine,atp] = synchronize(RecLine,SynPart,SymbolLen,GI,TimeSlot*(SymbolNum+1)+1);
    disp(['帧头位置：',num2str(atp)]);
    FrameRec_GI = reshape(FrameLine,SymbolLen,(SymbolNum+1)*TimeSlot);
    FrameRecGT = FrameRec_GI.';
    FrameRec = FrameRecGT(:,GI+1:end);                                     % 去CP
    Pilot = FrameRec(1:SymbolNum+1:end,:);                                 % 提导频
    
    if RLS==1
        % ----------------------------------------------   RLS均衡   ----------------------------------------------
        PilotAve = mean(Pilot);
        PilotAve = PilotAve - mean(PilotAve);
        Gain = (max(PilotBlock)-min(PilotBlock))/(max(PilotAve)-min(PilotAve));
        PilotAve = PilotAve*Gain;
        A = [PilotBlock,PilotBlock,PilotBlock,PilotBlock,PilotBlock,PilotBlock];
        B = [PilotAve,PilotAve,PilotAve,PilotAve,PilotAve,PilotAve];
        [a,b,c] = RLSEqualizer(A,B,Order,sigma);
        FrameLine = FrameLine*Gain;
        Res = filter(a,1,FrameLine);
        FrameRec_GI = reshape(Res,IFFTLength+GI,(SymbolNum+1)*TimeSlot);
        FrameRecGT = FrameRec_GI.';
        FrameRec = FrameRecGT(:,GI+1:end);                                 % 去CP
        FrameRec(1:SymbolNum+1:end,:) = [];                                % 去导频
        TimeSlotCor = fft(FrameRec,IFFTLength,2);
    elseif RLS==2
        % ---------------------------------------------- LS信道估计 ----------------------------------------------
        FrameRec(1:SymbolNum+1:end,:) = [];                                % 去导频
        FrameRec = fft(FrameRec,IFFTLength,2);
        OriPilot= fft(PilotBlock);
        PilotAve = mean(Pilot);
        Gain = (max(PilotBlock)-min(PilotBlock))/(max(PilotAve)-min(PilotAve));
        PilotAve = PilotAve*Gain;
        
        H_LS = fft(PilotAve)./OriPilot;
        TimeSlotCor = zeros(TimeSlot*SymbolNum,IFFTLength);
        for floor = 1:TimeSlot*SymbolNum
            TimeSlotCor(floor,:) = FrameRec(floor,:)./H_LS;
        end
    elseif RLS ==0
        FrameRec(1:SymbolNum+1:end,:) = [];                                % 去导频
        FrameRec = fft(FrameRec,IFFTLength,2);
        TimeSlotCor = FrameRec;
    end
    
    % 信号提取
    ExtractMap = TimeSlotCor(:,2:1+SubCarrier);                  
    Constellation = [];
    % CASE 1
    if item==1
        % 子载波情况分布
        for sub=1:63
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % 用户EVM计算
        EVMuser(cyc,1) = mean(EVMsub(cyc,1:21));
        EVMuser(cyc,2) = mean(EVMsub(cyc,22:42));
        EVMuser(cyc,3) = mean(EVMsub(cyc,43:63));
        if status ==1
            EVMName = sprintf('%s\\UserAllEvmSub%d.txt',Folder,Num);
            save(EVMName, 'EVMsub', '-ASCII', '-double');
        end
    % CASE 23
    elseif item==2
        for sub=1:21
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % 用户EVM计算
        EVMuser(cyc,1) = mean(EVMsub(cyc,1:21));
        EVMName = sprintf('%s\\User1EvmSub%d.txt',Folder,Num);
        save(EVMName, 'EVMsub', '-ASCII', '-double');
        
    elseif item==3
        for sub=22:42
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % 用户EVM计算
        EVMuser(cyc,2) = mean(EVMsub(cyc,22:42));
        EVMName = sprintf('%s\\User2EvmSub%d.txt',Folder,Num);
        save(EVMName, 'EVMsub', '-ASCII', '-double');
    elseif item==4
        % 子载波情况分布
        for sub=43:63
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % 用户EVM计算
        EVMuser(cyc,3) = mean(EVMsub(cyc,43:63));
        if status == 1
            EVMName = sprintf('%s\\User3EvmSub%d.txt',Folder,Num);
            save(EVMName, 'EVMsub', '-ASCII', '-double');
        end
    end
end
% EVMsub(5,:) = [];
if item==1
    AveEVM1 = mean(EVMuser(:,1));
    AveEVM2 = mean(EVMuser(:,2));
    AveEVM3 = mean(EVMuser(:,3));
    AveEVM = [AveEVM1,AveEVM2,AveEVM3];
elseif item==2
    AveEVM = mean(EVMuser(:,1));
elseif item==3
    AveEVM = mean(EVMuser(:,2));
elseif item==4
    AveEVM = mean(EVMuser(:,3));
end

disp(' ');
disp(['当前实验情况:',num2str(ceil(item/3)),'  ','用户编号:',num2str(mod(item,3))]);
disp(['平均EVM: ',num2str(AveEVM)]);
scatterplot(Constellation);
title('4QAM星座图');

figure();
if Cycles==1
    plot(EVMsub);
    xlabel('Subcarries');
    ylabel('EVM(dB)');
    title('子载波EVM分布');
else
    plot(mean(EVMsub));
    xlabel('Subcarries');
    ylabel('EVM(dB)');
    title('平均子载波EVM分布');
end