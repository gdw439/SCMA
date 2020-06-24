% DataIn:    待同步数据
% SynPart:  同步序列
% SymLen: 符号长度
% CPLen:   循环前缀长度
% SymNum:符号数目

% OutData: 同步后的传输数据
% location:  同步前数据点位置
function [OutData, location] = synchronize(DataIn, SynPart, SymLen, CPLen, SymNum)
    DataIn = DataIn - mean(DataIn);
    LineLen = length(DataIn);
    Raw = floor(LineLen/SymLen) - 1;
    
    CorrSave = zeros(1,SymLen);
    for item = 1:SymLen
        temp = DataIn(item:item+Raw*SymLen-1);
        temp_mat = reshape(temp,SymLen,Raw);
        Corr = temp_mat(1:CPLen,:).*conj(temp_mat(SymLen-CPLen+1:SymLen,:));
        CorrSum = sum(Corr);
        CorrSave(item) = mean(CorrSum);
    end
    [~,StartPoint] = max(CorrSave);
    
    DataAlignment = DataIn(StartPoint:end);
    DataAlignment = DataAlignment - mean(DataAlignment);
    
    DataPart = DataAlignment(1:SymNum*SymLen);
    DataMat = reshape(DataPart,SymLen,SymNum);
    DataMat = DataMat.';
    tem2 = conj(SynPart);
    
    CorrSave1 = zeros(1,SymNum);
    for num = 1:SymNum
        CorrSave1(num) = sum(DataMat(num,:).*tem2);
    end
    [~,StaSym] = max(CorrSave1);
    StaPoi = StaSym*SymLen+1
    EndPoi = StaPoi+SymLen*SymNum-1
    
    OutData = DataAlignment(StaPoi:EndPoi);
    OutData = OutData - mean(OutData);
    location = StartPoint + StaPoi;
    
    % ---------------- 调试部分 --------------------
%     figure();
%     subplot(2,1,1);
%     plot(CorrSave);
%     xlabel('位移量');
%     ylabel('相关幅度');
%     title('OFDM符号对齐相关峰');
%     subplot(2,1,2);
%     plot(CorrSave1);
%     xlabel('符号数');
%     ylabel('相关幅度');
%     title('起始符号搜索');
end