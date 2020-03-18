function [CEN_N,SIG_N,range]=newCENTER2B (numS,US,LS,DIM,CEN_P,p,z)
% [c,v,r] = newCENTER2B (numS,US(k),LS(k),k,C(ind,k),p,z)
%numS======THE NUMBER OF SPILT
%US========UPPER bound of corresponding spilt
%LS========LOWER bound of corresponding spilt
%  global z
%  global p   
%  R=zeros(numS)
 if numS==2
    R(1,:)=[LS,CEN_P];
    R(2,:)=[CEN_P,US];
 end

 for i=1:numS
 
 Z=[] ;
 DS=[];
 [Z ]=find(z(:,DIM)>=R(i,1)& (z(:,DIM)<=R(i,2)));
 DS=z(Z,DIM);
 CEN_N(i)= mean(DS);
  Dis(i)=max(CEN_N(i)-R(i,1),R(i,2)-CEN_N(i));
 SIG_N(i)=Dis(i)/2;
 
 end
 
 range=zeros(numS,2);
 range=R;