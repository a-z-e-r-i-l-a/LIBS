clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength_range=[200,800];
wavelength_resolution=1;
a_fraction=0.01;     %remove peaks below a_fraction of the maximum peak

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



F = dir('*.csv');
select_path = uigetdir;
wavelengths=[wavelength_range(1):wavelength_resolution:wavelength_range(2)];
for class=[1:length(F)]

name=F(class).name;

fileID = fopen(name);
A=textscan(fileID,' %s %s','Delimiter',',','EmptyValue',0);
fclose(fileID);
    for lambda=[1:length(A{1,1})]
        spectrum(lambda,1)=str2double(A{1,1}{lambda,1});
        [d,idx] = min( abs( wavelengths-spectrum(lambda,1) ) );
        spectrum(lambda,1)=wavelengths(idx);
        spectrum(lambda,2)=str2double(A{1,2}{lambda,1});
    end
    spectrum(:,2)=spectrum(:,2)/max(spectrum(:,2));  %normalization of intensities
                %writing the main spectrum vector
    [d,idx]=find(wavelengths==spectrum(:,1));
    main_spectrum=(zeros(length(wavelengths),1))';
    main_spectrum(idx)=spectrum(:,2);
                %removing peaks below 0.1 of max
                [d2,idx2]=find(main_spectrum<a_fraction);
                main_spectrum(idx2)=0;
                
                  
        Data_set{class,1}=name; labels{class}=name;
        Data_set{class,2}=class;
        for i=[3:length(main_spectrum)+2]
        Data_set{class,i}=main_spectrum(i-2);
        X_train{class,i-2}=main_spectrum(i-2);
        end
clear idx;clear idx2;clear main_spectrum;clear spectrum;
end
%Labels
%cell2csv('Data_set.data',X_train,',','w')

for i=[1:601]
    a(i)=X_train{30,i};
end
find(a~=0)

