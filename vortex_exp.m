close all;clear;clc;
format long

% ---------------------------------------- ����̨  ----------------------------------------
        RLS = 2;                                                           % 0 --�޾���   | 1 --RLS����  | 2 --LS�����������
        item = 1;                                                          % ʵ�鲽��(1--case1; 2~4--case2)
        status = 2;                                                        % 0 --����״̬ | 1 --ʵ��״̬
        Fs = 50e6;                                                         % ����Ƶ��
        Cycles = 1;                                                        % �ظ�����
        TimeSlot = 30;                                                     % ��������
        
        sigma = 0.6;                                                      % RLS ����[0.8:2]
        Order =3;
        
        Cut = 0;                                                           % ���б�־λ
        factor = 13.5;                                                     % ������
        
        Folder = '�û�ʵ��';                                                % ���Ŀ¼
        OveSamT = 4;                                                       % ����������
        OffLineData = '�û�ʵ��\ExpData1.mat';                            % ��������Ŀ¼
        
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
SymbolNum = 24;                                                            % OFDM ������
IFFTLength = 128;
GI = PrefixRatio*IFFTLength;
SubCarrier = IFFTLength/2-1;                                               % ����������ز�
BitSymbol  = modebit(curmode);
EvmRate = zeros(1,length(Cycles));                                         % EVM
EVMsub = zeros(Cycles,SubCarrier);
EVMuser= zeros(Cycles,3);
SymbolLen = IFFTLength+GI;
TimeSlotLen = (SymbolNum+1)*SymbolLen;
FrameLon = TimeSlot*SubCarrier*SymbolNum*BitSymbol;                        % ��������
UserTotal = TimeSlot*SubCarrier*SymbolNum/3;
disp(['Total simulation: ',num2str(FrameLon)]);
disp(' ');

% ͬ������
N = 63;
n = 1:1:N;
nzc = exp( -1i*pi*34*(n-1).*n./N);                                         % ZC����
Nzc = [nzc,nzc(1)];
NzcJ = conj(Nzc);
SynSeq = zeros(1,IFFTLength);
SynSeq(2:1+length(Nzc)) = Nzc;
SynSeq(end:-1:end-length(Nzc)+1) = NzcJ(1:end);                            % �����
SynBlock = ifft(SynSeq);
SynPart = [SynBlock(end-GI+1:end),SynBlock];                               % ���CP

% ��Ƶ
N = SubCarrier-1;
n = 1:1:N;
nzc = exp( -1i*pi*25*(n-1).*n./N);
Nzc = [nzc,nzc(1)];
PilotSeq = zeros(1,IFFTLength);
PilotSeq(2:1+length(Nzc)) = Nzc;
NzcJ = conj(Nzc);
PilotSeq(end:-1:end-length(Nzc)+1) = NzcJ(1:end);   % ����
PilotBlock = ifft(PilotSeq);
% -------------------------------------------------------  ���Ͷ�  ------------------------------------------------------- 

if exist('�û�ʵ��','dir') == 0
   mkdir('�û�ʵ��'); 
end

if exist('�û�ʵ��\Info.txt','file') == 0  
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
        
        % ʱ��任
        FrameBox = zeros(TimeSlot*SymbolNum,IFFTLength);
        FrameBox(:,2:SubCarrier+1) = FrameMapT;
        FrameBox(:,end:-1:end-SubCarrier+1) = conj(FrameMapT);
        FrameBox = ifft(FrameBox,IFFTLength,2);
        FrameTran = zeros(TimeSlot*(SymbolNum+1),IFFTLength);
            
        Peak = max(max(abs(FrameBox)));                                    % ��һ����ֵ������ͻ��ë��
        SynPart = Peak*SynPart./max(SynPart)./3;
        PilotBlock = Peak*PilotBlock./max(PilotBlock)./3;
                
        % ��Ƶ��ǰ׺
        for slot = 0:TimeSlot-1
            FrameTran(slot*(SymbolNum+1)+1,:) = PilotBlock;
            FrameTran((slot*(SymbolNum+1)+2):(slot*(SymbolNum+1)+1+SymbolNum),:) = FrameBox(slot*SymbolNum+1:(slot+1)*SymbolNum,:);
        end
        FramePort = [FrameTran(:,IFFTLength-GI+1:end),FrameTran];
        FramePort = reshape(FramePort.',1,TimeSlotLen*TimeSlot);

        % ͬ������
        TranLine = [SynPart,FramePort];               
        SignalBackUp = TranLine;

        % ������[4 times]
        Len = length(TranLine);
        TranLine = oversamp(TranLine, Fs, Fs*OveSamT, 0, 0,(OveSamT - 1)*Len);
        TranLine = real(TranLine);

        % �źż���  
        if Cut == 1
            factor = 10^(factor/10);
            signalpower = mean(TranLine.^2);
            ClippingThreshold = sqrt(factor*signalpower);                  % ���������
            TranLine(TranLine> ClippingThreshold) = ClippingThreshold;
            TranLine(TranLine<-ClippingThreshold) =-ClippingThreshold;
        end
                
        % �洢�����ź�
        TranLine = TranLine.';
        DataName = sprintf('%s\\User%d Signal.txt',Folder,user);
        save(DataName, 'TranLine', '-ASCII', '-double');
    end
    
    % �洢����ͼ
    FrameMapr = real(FrameMap);
    FrameMapi = imag(FrameMap);
    MapName = sprintf('%s\\FrameMapr.txt',Folder);
    save(MapName, 'FrameMapr', '-ASCII', '-double');
    MapName = sprintf('%s\\FrameMapi.txt',Folder);
    save(MapName, 'FrameMapi', '-ASCII', '-double');
    DataLen = length(TranLine);
    disp(['����������Ϊ��',num2str(DataLen)]);
end


% ---------------------------------------- ������������  ----------------------------------------
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
    disp('----------------------��ǰΪ����״̬----------------------')
    DataName = sprintf('%s\\User1 Signal.txt',Folder);
    UserSignal1 = importdata(DataName);
    
    DataName = sprintf('%s\\User2 Signal.txt',Folder);
    UserSignal2 = importdata(DataName);
    
    DataName = sprintf('%s\\User3 Signal.txt',Folder);
    UserSignal3 = importdata(DataName);
    
    disp('���ݵ������!');
    TranLine = UserSignal1 + UserSignal2 + UserSignal3;
    TranLine = TranLine.';
elseif status == 1                                                         % ���ݵ���
    disp('----------------------��ǰΪʵ��״̬----------------------')
    Point = 2e6;
    OSC = OSCInit();        % ����ʾ����
    OSCPoint(OSC,Point);    % ���ò�����
elseif status == 2
    disp('----------------------��ǰΪ����״̬----------------------')
end

for cyc = 1:Cycles
    disp(['��ǰѭ��:',num2str(cyc)]);
    if status == 1                                                         % ʵ��״̬
        TranLine = OSCRead(OSC);
        disp(['���ȣ�',num2str(length(TranLine))]);
        
        DataName = sprintf('%s\\Info.txt',Folder);
        Num = load(DataName);
        Num = Num + 1;
        filename = sprintf('%s\\ExpData%d',Folder,Num);
        save(strcat(filename,'.mat'),'TranLine');
        save(DataName, 'Num', '-ASCII', '-double');
        disp('�ռ����ݴ洢��ϣ�');
    elseif status == 2                                                     % ����״̬
%         OffLineData = sprintf('�û�ʵ��\\ExpData%d.mat',cyc);
        as = load(OffLineData);
        TranLine = as.TranLine;
    end

    % ----------------------------------------------      LPF     ----------------------------------------------
    lpf = zeros(1,length(TranLine));
    lpf(1:1:ceil(end/2/OveSamT)) = 1;
    lpf(end:-1:ceil(end-end/2/OveSamT+1)) = 1;
    TranLine = real(ifft(lpf.*fft(TranLine)));
    
%     TranLine = LMS(TranLine);
    % ������
    RecLine = TranLine(1:OveSamT:end);
    
    RecLine = RecLine - mean(RecLine);
    [FrameLine,atp] = synchronize(RecLine,SynPart,SymbolLen,GI,TimeSlot*(SymbolNum+1)+1);
    disp(['֡ͷλ�ã�',num2str(atp)]);
    FrameRec_GI = reshape(FrameLine,SymbolLen,(SymbolNum+1)*TimeSlot);
    FrameRecGT = FrameRec_GI.';
    FrameRec = FrameRecGT(:,GI+1:end);                                     % ȥCP
    Pilot = FrameRec(1:SymbolNum+1:end,:);                                 % �ᵼƵ
    
    if RLS==1
        % ----------------------------------------------   RLS����   ----------------------------------------------
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
        FrameRec = FrameRecGT(:,GI+1:end);                                 % ȥCP
        FrameRec(1:SymbolNum+1:end,:) = [];                                % ȥ��Ƶ
        TimeSlotCor = fft(FrameRec,IFFTLength,2);
    elseif RLS==2
        % ---------------------------------------------- LS�ŵ����� ----------------------------------------------
        FrameRec(1:SymbolNum+1:end,:) = [];                                % ȥ��Ƶ
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
        FrameRec(1:SymbolNum+1:end,:) = [];                                % ȥ��Ƶ
        FrameRec = fft(FrameRec,IFFTLength,2);
        TimeSlotCor = FrameRec;
    end
    
    % �ź���ȡ
    ExtractMap = TimeSlotCor(:,2:1+SubCarrier);                  
    Constellation = [];
    % CASE 1
    if item==1
        % ���ز�����ֲ�
        for sub=1:63
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % �û�EVM����
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
        
        % �û�EVM����
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
        
        % �û�EVM����
        EVMuser(cyc,2) = mean(EVMsub(cyc,22:42));
        EVMName = sprintf('%s\\User2EvmSub%d.txt',Folder,Num);
        save(EVMName, 'EVMsub', '-ASCII', '-double');
    elseif item==4
        % ���ز�����ֲ�
        for sub=43:63
            Linew = ExtractMap(:,sub);
            Linew = Linew.';
            MeanPower = Linew.*conj(Linew);
            SignalMean = sqrt(mean(MeanPower));
            Linew = Linew./SignalMean;
            EVMsub(cyc,sub) = evm_cal(FrameMapT(:,sub).',Linew);
            Constellation = [Constellation,Linew];
        end
        
        % �û�EVM����
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
disp(['��ǰʵ�����:',num2str(ceil(item/3)),'  ','�û����:',num2str(mod(item,3))]);
disp(['ƽ��EVM: ',num2str(AveEVM)]);
scatterplot(Constellation);
title('4QAM����ͼ');

figure();
if Cycles==1
    plot(EVMsub);
    xlabel('Subcarries');
    ylabel('EVM(dB)');
    title('���ز�EVM�ֲ�');
else
    plot(mean(EVMsub));
    xlabel('Subcarries');
    ylabel('EVM(dB)');
    title('ƽ�����ز�EVM�ֲ�');
end