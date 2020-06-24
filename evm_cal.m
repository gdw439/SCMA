% ����EVMֵ
% valueR: ���ն��ź�;
% valueT: ���Ͷ��ź�;

function evm = evm_cal(valueT,valueR)
    if length(valueR)~=length(valueT)
        error('���ݳ��Ȳ�һ�£�');
    end
    Verr = valueT-valueR;
    Eerr =mean(Verr.*conj(Verr));
    Estd = mean(abs(valueT).^2);
    evm = Eerr/Estd;
    evm = 10*log10(evm);
end