%% load the merkmal
load('can_vmps_yawrate_chassis.mat', 'Effwert')
Daten(:,1)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'Effwert')
Daten(:,2)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'Effwert')
Daten(:,3)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'Effwert')
Daten(:,4)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'Effwert')
Daten(:,5)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'Effwert')
Daten(:,6)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'Effwert')
Daten(:,7)=Effwert(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'Effwert')
Daten(:,8)=Effwert(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'Effwert')
Daten(:,9)=Effwert(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'Effwert')
Daten(:,10)=Effwert(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'FFTwert_max')
Daten(:,11)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'FFTwert_max')
Daten(:,12)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'FFTwert_max')
Daten(:,13)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'FFTwert_max')
Daten(:,14)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'FFTwert_max')
Daten(:,15)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'FFTwert_max')
Daten(:,16)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'FFTwert_max')
Daten(:,17)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'FFTwert_max')
Daten(:,18)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'FFTwert_max')
Daten(:,19)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'FFTwert_max')
Daten(:,20)=FFTwert_max(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'FFTwert_mean')
Daten(:,21)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'FFTwert_mean')
Daten(:,22)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'FFTwert_mean')
Daten(:,23)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'FFTwert_mean')
Daten(:,24)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'FFTwert_mean')
Daten(:,25)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'FFTwert_mean')
Daten(:,26)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'FFTwert_mean')
Daten(:,27)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'FFTwert_mean')
Daten(:,28)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'FFTwert_mean')
Daten(:,29)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'FFTwert_mean')
Daten(:,30)=FFTwert_mean(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'FFTwert_varianz')
Daten(:,31)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'FFTwert_varianz')
Daten(:,32)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'FFTwert_varianz')
Daten(:,33)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'FFTwert_varianz')
Daten(:,34)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'FFTwert_varianz')
Daten(:,35)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'FFTwert_varianz')
Daten(:,36)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'FFTwert_varianz')
Daten(:,37)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'FFTwert_varianz')
Daten(:,38)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'FFTwert_mean')
Daten(:,39)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'FFTwert_mean')
Daten(:,40)=FFTwert_varianz(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'Woelbung')
Daten(:,41)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'Woelbung')
Daten(:,42)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'Woelbung')
Daten(:,43)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'Woelbung')
Daten(:,44)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'Woelbung')
Daten(:,45)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'Woelbung')
Daten(:,46)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'Woelbung')
Daten(:,47)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'Woelbung')
Daten(:,48)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'Woelbung')
Daten(:,49)=Woelbung(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'Woelbung')
Daten(:,50)=Woelbung(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'schiefe')
Daten(:,51)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'schiefe')
Daten(:,52)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'schiefe')
Daten(:,53)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'schiefe')
Daten(:,54)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'schiefe')
Daten(:,55)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'schiefe')
Daten(:,56)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'schiefe')
Daten(:,57)=schiefe(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'schiefe')
Daten(:,58)=schiefe(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'schiefe')
Daten(:,59)=schiefe(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'schiefe')
Daten(:,60)=schiefe(:,1);%Eff(luft strom koeperschall)
%% merkmale
load('can_vmps_yawrate_chassis.mat', 'varianz')
Daten(:,61)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_vmps_rollrate_chassis.mat', 'varianz')
Daten(:,62)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_vmps_pitchrate_chassis.mat', 'varianz')
Daten(:,63)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm73_rollrate_chassis.mat', 'varianz')
Daten(:,64)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm73_acceleration_chassis_zneg.mat', 'varianz')
Daten(:,65)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm72_pitchrate_chassis.mat', 'varianz')
Daten(:,66)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm72_acceleration_chassis_zneg.mat', 'varianz')
Daten(:,67)=varianz(:,1);%Eff(luft strom koeperschall)
load('can_mm71_yawrate_chassis.mat', 'varianz')
Daten(:,68)=varianz(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_rollnoise_right.mat', 'varianz')
Daten(:,69)=varianz(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'varianz')
Daten(:,70)=varianz(:,1);%Eff(luft strom koeperschall)
load('analog_daqmx_acceleration_supension_z_right.mat', 'Label')
Daten(:,71)=Label;


%% normal
for d=1:70
Daten(:,d)=(Daten(:,d)-mean(Daten(:,d)))/std(Daten(:,d));
end
%%
%n=3;
%klassenvector=[79 9];
%Daten=sortrows(Daten,51);
%merkmalvek=addonalgo(Daten(:,1:50),n,klassenvector);

 %traindata=D(:,merkmalvek);
 %testdata=Daten(:,merkmalvek);
%%
%Target=Testdata(:,13);%EKR_Bayes=1-sum(Target~=yfit)/length(Target)
%yfit = trainedClassifier.predictFcn(testdata(:,1:3));
