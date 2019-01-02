%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in  this model amygdala have output from thalamous  
%in this method we train a network with conditional stimuli
%then evaluation of function of system estimate and Emotional asignal prepare
%this estimate is in manner that responce of system to new stimuli was supposed is correct and the responce was supposed as a valid target 
%agai the system costruct with new stimula and target 
%for this system output for conditional stimula was computed and the output contrast with target and emotional signl prepar 
%insted use of whole of train or condition stimula we use 1 of them   

%%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDDDDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test size 1800 t0 2000
for sim=1:10
UP2=311+sim-1%%%yaer 1965
LP2=308
TSS= yearsnn_new(LP2:UP2,:);
in=2
MM=[];

for j=1+in:UP2-LP2+1
    MM(j-in,:) = [ TSS(j-2,2) TSS(j-1,2) TSS(j,2)];
end
aa1=[];
aa1=MM;
U=[aa1];
aan=aa1;
N_TEST=size(MM,1);
p=in;
ST=aan(:,1:p);
Atht=max(ST')';
SST=[ST,Atht];
PT=SST;
PPT=[ST];



Ea=evalfis(PPT,fismatA);





EET(sim)=Ea(end);
% EET1 =ET.*(max(t)-min(t))+min(t) ;

yearsnn_new(UP2,2)=EET(sim);
yearsnn_new(UP2+1,2)=EET(sim);
end