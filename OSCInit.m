function OSC = OSCInit()
    visaResourceString = 'USB0::2391::6032::MY54390336::0::INSTR';
    obj1 = instrfind('Type', 'visa-usb', 'RsrcName', visaResourceString, 'Tag', '');
    if isempty(obj1)
        obj1 = visa('AGILENT', visaResourceString);
    else
        fclose(obj1);
        obj1 = obj1(1);
    end
    set(obj1, 'InputBufferSize', 8e6);
    set(obj1, 'OutputBufferSize', 512);
    fopen(obj1);
    OSC = obj1;
end