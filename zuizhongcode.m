clc
clear

%% 提取数据信息
info_3 = hdf5info('../m1_88-200.hdf5');     % label的数据
%load('../../conv_tdms/m1_88-200.mat');
L=length(info_3.GroupHierarchy.Groups);     % label的个数或者长度
fre_data=51200;                             % 51200采样频率

for i=1:L                                   % 计算总的label
    label  = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(1));  % 读取第i组label
    time_3 = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(2));  % 读取对应时间
    %name_3 = info_3.GroupHierarchy.Groups(i).Name;
    if  ~all(label==0)                         %判断是否label存在错误情况
        index=find(label>0);
        sumlabel(index)=label(index);
        %ss=time_3(index);
        %eval(['time_value',num2str(i),'=','ss',';'])
          
    end
end
        sumlabel(1,length(label))=0;         %最终的label
        xx=time_3(1):(1/51200):time_3(end);  %插值label
        sumlabel=interp1(time_3,sumlabel,xx, 'nearest'); %插值后的label
       
%% 设置window窗
        window=0.02;             %20ms
        slice_data=fre_data*window;  %window 长度
        N=round((time_3(end)-time_3(1))/window);  % N is the number of slices per window sencond
        Label(N,1)=0;
        Target=Label;
        sumlabel(1,N*slice_data)=0;  %插值后的label
%sumlabel=sumlabel';
for ii=1:10
    Data=0;      %初始化
    str1=sprintf('data%d',ii);
    str2=sprintf( '.mat');
    str=[ str1 str2 ];
    load(sprintf(str));
%% analog data import
    Data=data.data;% original signal
    time=data.time;% original time
    %fre_data=round(length(Data)/(time(end)-time(1)));


%% 对数据进行interpolation
    yy=time(1):(1/51200):time(end);  %插值的时间
    Data=interp1(time,Data,yy, 'spline');
    Data(1,N*slice_data)=0;   %插值后的数据
    %Data=Data';


%% 计算merkmal


%time_3 = hdf5read(info_3.GroupHierarchy.Groups(1).Datasets(2)); % time of label
%fre_label=round(length(label)/(time_3(end)-time_3(1))); % frequency of label
        fre_label=51200;
        slice_label=fre_label*window;  % index of the label per window
%N=round(time_3(end)-time_3(1));
        for iii=1:N                                        % N slices
            index_data=(iii-1)*slice_data+1:iii*slice_data;
            aa=Data(index_data);
            index_label=(iii-1)*slice_label+1:iii*slice_label;
            
           if ~all(sumlabel(index_label)==0)
                %eval(['data_bad',num2str(ii),'=','aa',';']);
                Effwert(iii,1:2)=[Eff(aa),1];
                varianz(iii,1:2)=[var(aa),1];
                schiefe(iii,1:2)=[skewness(aa),1];
                %statismerkmale(ii,4)=mean(aa);
               [Fre,Amp] = frequencySpectrum(aa,fre_data);
                FFTwert_mean(iii,1:2)=[mean(Amp(2:end)),1];
                FFTwert_varianz(iii,1:2)=[var(Amp(2:end)),1];
                FFTwert_max(iii,1:2)=[max(Amp(2:end)),1];
                %statismerkmale(ii,5)=kurtosis(aa);
                %statismerkmale(ii,6)=max(aa);
                %Varianz(ii,3)=var(aa);
                %Schiefe(ii,3)=skewness(aa);
                Woelbung(iii,1:2)=[kurtosis(aa),1];
                Target(iii,1)=1;
           else 
                %eval(['data_gut',num2str(ii),'=','aa',';']);
                Effwert(iii,1:2)=[Eff(aa),0];
                varianz(iii,1:2)=[var(aa),0];
                schiefe(iii,1:2)=[skewness(aa),0];
                [Fre,Amp] = frequencySpectrum(aa,fre_data);
                FFTwert_mean(iii,1:2)=[mean(Amp(2:end)),0];
                FFTwert_varianz(iii,1:2)=[var(Amp(2:end)),0];
                FFTwert_max(iii,1:2)=[max(Amp(2:end)),0];
                %statismerkmale(ii,4)=mean(aa);
                %[Fre,Amp] = frequencySpectrum(aa,200);
                % FFTwert(ii,1)=Amp(48);
                %statismerkmale(ii,5)=kurtosis(aa);
                %statismerkmale(ii,6)=max(aa);
                %Varianz(ii,3)=var(aa);
                %Schiefe(ii,3)=skewness(aa);
                Woelbung(iii,1:2)=[kurtosis(aa),0];
                Target(iii,1)=0;
          
   % figure()
    %plot(time(index(1):index(end)),ss);
    %hold on
    %plot(time_3(index_3(1)-100:index_3(end)+100),data_3(index_3(1)-100:index_3(end)+100))
    %legend('analog','label')
           end 
           
        end
          % save Varianz Varianz  
%% merkmalraum darstellung

 Daten(:,ii*7-6:ii*7)=[Effwert(:,1),varianz(:,1),schiefe(:,1),FFTwert_mean(:,1),FFTwert_varianz(:,1),FFTwert_max(:,1),Woelbung(:,1)];
end
%save analog_suspension_z_right Merkmal
%save analog_daqmx_acceleration_supension_z_right
