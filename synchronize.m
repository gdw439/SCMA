% DataIn:    ��ͬ������
% SynPart:  ͬ������
% SymLen: ���ų���
% CPLen:   ѭ��ǰ׺����
% SymNum:������Ŀ

% OutData: ͬ����Ĵ�������
% location:  ͬ��ǰ���ݵ�λ��
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
    
    % ---------------- ���Բ��� --------------------
%     figure();
%     subplot(2,1,1);
%     plot(CorrSave);
%     xlabel('λ����');
%     ylabel('��ط���');
%     title('OFDM���Ŷ�����ط�');
%     subplot(2,1,2);
%     plot(CorrSave1);
%     xlabel('������');
%     ylabel('��ط���');
%     title('��ʼ��������');
end