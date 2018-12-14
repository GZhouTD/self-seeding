clc;clear
close all

ssnum = 1:100;
snum = 51;
R56 = 100;
filter = 3;
foldname =  ['genesis/dump_1_XLBEG_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter)];
data = cell(length(ssnum),1);
for ss = 1:length(ssnum)
    ss
    [data1,info1] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis.out']);
    if ssnum(ss)<=50
        [data2,info2] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis2_15.out']);
    else
        [data2,info2] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis3_15.out']);
    end
    power1 = data1.power(end,:)/1e9;
    power2 = data2.power(end,:)/1e9;
    t1 = data1.t/3e8/1e-15;
    t2 = data2.t/3e8/1e-15;
    ff1 = data1.freq;
    ff2 = data2.freq;
    freq1=299792458/info1.lambda*6.62606957e-34/1.60217657e-19;
    freq2=299792458/info2.lambda*6.62606957e-34/1.60217657e-19;
    df1 = (ff1-freq1)/freq1;
    df2 = (ff2-freq2)/freq2;
    spec1 = data1.spectrum(end,:);
    spec2 = data2.spectrum(end,:);
    result.t1 = t1;
    result.t2 = t2;
    result.ff1 = ff1;
    result.ff2 = ff2;
    result.df1 = df1;
    result.df2 = df2;
    result.freq1 = freq1;
    result.freq2 = freq2;
    result.power1 = power1;
    result.power2 = power2;
    result.spec1 = spec1;
    result.spec2 = spec2;
    data{ss,1} = result;
end
save('two_color_genesis2_15.mat','data','-v7.3')