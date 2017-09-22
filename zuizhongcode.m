clc
clear

%% ��ȡ������Ϣ
info_3 = hdf5info('../m1_88-200.hdf5');     % label������
%load('../../conv_tdms/m1_88-200.mat');
L=length(info_3.GroupHierarchy.Groups);     % label�ĸ������߳���
fre_data=51200;                             % 51200����Ƶ��

for i=1:L                                   % �����ܵ�label
    label  = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(1));  % ��ȡ��i��label
    time_3 = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(2));  % ��ȡ��Ӧʱ��
    %name_3 = info_3.GroupHierarchy.Groups(i).Name;
    if  ~all(label==0)                         %�ж��Ƿ�label���ڴ������
        index=find(label>0);
        sumlabel(index)=label(index);
        %ss=time_3(index);
        %eval(['time_value',num2str(i),'=','ss',';'])
          
    end
end
        sumlabel(1,length(label))=0;         %���յ�label
        xx=time_3(1):(1/51200):time_3(end);  %��ֵlabel
        sumlabel=interp1(time_3,sumlabel,xx, 'nearest'); %��ֵ���label
       
%% ����window��
        window=0.02;             %20ms
        slice_data=fre_data*window;  %window ����
        N=round((time_3(end)-time_3(1))/window);  % N is the number of slices per window sencond
        Label(N,1)=0;
        Target=Label;
        sumlabel(1,N*slice_data)=0;  %��ֵ���label
%sumlabel=sumlabel';
for ii=1:10
    Data=0;      %��ʼ��
    str1=sprintf('data%d',ii);
    str2=sprintf( '.mat');
    str=[ str1 str2 ];
    load(sprintf(str));
%% analog data import
    Data=data.data;% original signal
    time=data.time;% original time
    %fre_data=round(length(Data)/(time(end)-time(1)));


%% �����ݽ���interpolation
    yy=time(1):(1/51200):time(end);  %��ֵ��ʱ��
    Data=interp1(time,Data,yy, 'spline');
    Data(1,N*slice_data)=0;   %��ֵ�������
    %Data=Data';


%% ����merkmal


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
