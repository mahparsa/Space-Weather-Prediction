% What(k,:)=inv(X'*PHIE(k,:)'*PHIE(k,:)*X)*X'*PHIE(k,:)'*PHIE(k,:)*yn
function [LOSS,W,MSEE]= calLOSSG6(CEN,VAR,R);
global z;
global yn;
global p;
global N;
global x;
a=0.5;
alpha=0.000000001;
xx=zeros(N,p+1);
I=p;
What=[];
TEMP=[];
MEIYU=[];
PHIE=[];
X=[];
% SIG=(1/sqrt(2))./VAR 
SIG=1./(sqrt(2)*VAR);
for l=1:R
        
    for k=1:I
            
        TEMP(:,k)= (z(:,k)-CEN(l,k))*SIG(l,k);

    end    
    TEMP=TEMP.*TEMP ; 
    MEIYU(l,:)=sum(TEMP,2)'; 

end 
Z=sum(MEIYU,1);

for l=1:R
     PHIE(l,:)= MEIYU(l,:) ./ Z;
end

for kk=1:R
     
      for k=1:p+1
          
          if k==1  
          
              xx(:,k)=PHIE(kk,:);
        
          else
          
              xx(:,k) = x(:,k-1).*PHIE(kk,:)';
          
          end
       
      end
 X=[X,xx]; 
  end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for kk=1:M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GLOBALLY WEIGHTED squre
What=inv( X'*X+alpha* eye(R*(p+1),R*(p+1)))*X'*yn;
yhat =X * What;
Y=yhat;
E=yn - Y;
% [Imin,Imax,PMAX,PMIN]=PEAKMINMAX(N,yn);
% for i=1:size(Imax)
% EMAX(i)=yn(i)-Y(i);
% end
% 
% for i=1:size(Imin)
% EMIN(i)=yn(i)-Y(i);
% end
for kk=1:R
%     loss(kk)=sum(mse(EMIN)*mse(EMAX)*PHIE(kk,:)'.* E.^2)%RMSE=0.116
%     loss(kk)=mse(PHIE(kk,:))*mse(E) %%%%RMSE=0.1253
%     loss(kk)=sum(PHIE(kk,:)'.*E.^2)%%%%RMSE=0125
%     loss(kk)=mse(E)+mse(EMAX)+mse(EMIN)RMSE=0.1008
%     loss(kk)=sum((1./PHIE(kk,:))'* sum(E.^2))%%%%RMSE=0.1631
%     loss(kk)=sum((1./PHIE(1,:))'.*E.^2)%%% RMSE=0.1247
%     loss(kk)=mse((1./PHIE(1,:))'.*E.^2)%%% RMSE=0.1247
%       loss(kk)=mse((1./PHIE(1,:))'.*E) %%%RMSE=0.0282
      loss(kk)=mse(PHIE(1,:)'.*(1./(E.^2)));%%%RMSE=0.016
%       loss(kk)=mse(PHIE(1,:)'.*(1./(E.^2)))/(mse(EMIN)+mse(EMAX));%%%%RMSE=0.0098

end 
 LOSS=loss;   
 W=What;
%  MSEE=mse(E)%%%RMSE=0.1
%  MSEE=mse(EMIN)+mse(EMAX);%%%RMSE=0.1