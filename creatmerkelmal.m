clear ;
%% 提取数据信息
info_3 = hdf5info('../m1_88-200.hdf5');
load('../../conv_tdms/m1_88-200.mat');
%% analog data import
L=length(info_3.GroupHierarchy.Groups);%数据的个数或者长度
data=can_vmps_acceleration_chassis_z.data;% original signal
time=can_vmps_acceleration_chassis_z.time;% original time
fre_data=round(length(data)/(time(end)-time(1)));
window=0.02;
slice_data=round(fre_data*window);  %window 长度
N=round((time(end)-time(1))/window);  % N is the number of slices per window sencond
Label(N,1)=0;
Target=Label;
sumlabel(:)=0;
data(slice_data*N+1,1)=0;
for i=1:L
    label = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(1));
    %time_3 = hdf5read(info_3.GroupHierarchy.Groups(i).Datasets(2));
    %name_3 = info_3.GroupHierarchy.Groups(i).Name;
    
    if  ~all(label==0)
        index=find(label>0);
        sumlabel(index)=label(index);
        %ss=time_3(index);
        %eval(['time_value',num2str(i),'=','ss',';'])
          
    end
end
sumlabel(1,length(label))=0;
sumlabel=sumlabel';

time_3 = hdf5read(info_3.GroupHierarchy.Groups(1).Datasets(2)); % time of label
fre_label=round(length(label)/(time_3(end)-time_3(1))); % frequency of label
slice_label=round(fre_label*window);  % index of the label per window
%N=round(time_3(end)-time_3(1));
        for ii=1:N                                        % N slices
            index_data=(ii-1)*slice_data+1:ii*slice_data;
            aa=data(index_data);
            index_label=(ii-1)*slice_label+1:ii*slice_label;
            
           if ~all(sumlabel(index_label)==0)
                %eval(['data_bad',num2str(ii),'=','aa',';']);
                Effwert(ii,1:2)=[Eff(aa),1];
                varianz(ii,1:2)=[var(aa),1];
                schiefe(ii,1:2)=[skewness(aa),1];
                %statismerkmale(ii,4)=mean(aa);
               [Fre,Amp] = frequencySpectrum(aa,fre_data);
                FFTwert_mean(ii,1:2)=[mean(Amp(2:end)),1];
                FFTwert_varianz(ii,1:2)=[var(Amp(2:end)),1];
                FFTwert_max(ii,1:2)=[max(Amp(2:end)),1];
                %statismerkmale(ii,5)=kurtosis(aa);
                %statismerkmale(ii,6)=max(aa);
                %Varianz(ii,3)=var(aa);
                %Schiefe(ii,3)=skewness(aa);
                Woelbung(ii,1:2)=[kurtosis(aa),1];
                Target(ii,1)=1;
           else 
                %eval(['data_gut',num2str(ii),'=','aa',';']);
                Effwert(ii,1:2)=[Eff(aa),0];
                varianz(ii,1:2)=[var(aa),0];
                schiefe(ii,1:2)=[skewness(aa),0];
                [Fre,Amp] = frequencySpectrum(aa,fre_data);
                FFTwert_mean(ii,1:2)=[mean(Amp(2:end)),0];
                FFTwert_varianz(ii,1:2)=[var(Amp(2:end)),0];
                FFTwert_max(ii,1:2)=[max(Amp(2:end)),0];
                %statismerkmale(ii,4)=mean(aa);
                %[Fre,Amp] = frequencySpectrum(aa,200);
                % FFTwert(ii,1)=Amp(48);
                %statismerkmale(ii,5)=kurtosis(aa);
                %statismerkmale(ii,6)=max(aa);
                %Varianz(ii,3)=var(aa);
                %Schiefe(ii,3)=skewness(aa);
                Woelbung(ii,1:2)=[kurtosis(aa),0];
                Target(ii,1)=0;
          
   % figure()
    %plot(time(index(1):index(end)),ss);
    %hold on
    %plot(time_3(index_3(1)-100:index_3(end)+100),data_3(index_3(1)-100:index_3(end)+100))
    %legend('analog','label')
           end 
           
           Label=Target|Label;
        end
          % save Varianz Varianz  
%% merkmalraum darstellung
 Merkmal(:,1)=Effwert(:,1);
 Merkmal(:,2)=varianz(:,1);
 Merkmal(:,3)=schiefe(:,1);
 Merkmal(:,4)=FFTwert_mean(:,1);
 Merkmal(:,5)=FFTwert_varianz(:,1);
 Merkmal(:,6)=FFTwert_max(:,1);
 Merkmal(:,7)=Woelbung(:,1);
 Label=Label+0;
%save can_vmps_acceleration_z Merkmal
%save can_vmps_acceleration_chassis_z
