%%%NORMALIZE IN [0 1]scale
% IN THIS MTHOD DO NOT USE CLUSTERING AT FIRST
% N =220  %%N show the sample data(training data)
% Nt=29%%Nt show the sample data(test  data)
global z;
global x;
% global xx
global lossORD;
global N;
global p;
global yn;


TRAIN_ONE=[];
for i=145:168
    TRAIN_ONE=[TRAIN_ONE,TRAIN1(i,:) ];
end
TEST_ONE=[];
for ii=1:24
    TEST_ONE=[TEST_ONE,TEST(ii,:) ];
end
UP1= size(TRAIN_ONE,2);
LP1=1;
LP2=1;
UP2=size(TEST_ONE,2);
TS=[];
TSS=[];
TS = TRAIN_ONE';
TSS= TEST_ONE';
M=[];
MM=[];
h=5;
in=2;


for i=in*h+1:UP1-LP1+1
    A(i-in*h,:) = [ TS(i-2*h,1) TS(i-h,1) TS(i,1)];
end
for j=1+in*h:UP2-LP2+1
    AA(j-in*h,:) = [   TSS(j-2*h,1) TSS(j-h,1) TSS(j,1)];
end



DATA=[A;AA];
Nt=size(DATA,1);
p=size(DATA,2)-1;
U = zeros(Nt,1:p)
aa1=DATA(:,1:p)
U=[aa1]
minU=min(U)
maxU=max(U)
for i=1:size(U,2)
    aan(:,i)=(U(:,i)- minU(i))./(maxU(i)-minU(i))
end

p=2
N=size(A,1);
NT=size(AA,1);
z=[];
x=[];
yn=[];
y=[];
z=aan(1:N,1:p);
x=aan(1:N,1:p);
y=DATA(1:N,p+1)
yn=(y- min(y))./( max(y)-min(y))


% yn=aan(1:N,p+1);
% y=A(:,p+1);
Sratio = 0.3; %% the splitting ratio that determine each partition devide to howmany partition 
alpha=0.00000001;
varR=0.5;
xx=zeros(N,p+1);
perf=[];
TItteration=100;
perf =ones(1,TItteration);
Itteration=1;
perf(Itteration)=1;
C=[];
V=[];
R=[];
X=[];
tic;
while ((perf(Itteration)>=0.0000000001)&&(Itteration<TItteration))

if Itteration==1
   C=mean(z);
   V=(max(z)-min(z));
   d=V-C;
   V=d./(sqrt(2)*sqrt(-log(0.4)));
   M=1;
else
    
    M=Law;
   
end % %%end related to Itteration==1
%%%%%%%%%%%%C and V are the current CENTER and SPEARD 
[LOSS,W]= calLOSSG6(C,V,M);
lossORD=LOSS;
What=W;
[ con, ind ] = min(lossORD);
CD=cell(1,p);
VD=cell(1,p);
% R=cell{1,p}%%%%%%EACH CELL USE FOR KEEPING THE RANGE OF EACH CENTER FOR this DIMENSION

spilt=2*ones(p,1);
numS=2;
if M==1
    for l=1:p
    R(l,:,M)=[min(z(:,l));max(z(:,l))];%%%%%%%%%%%%%%%%%%%%%range each center
    end
end
   
for k=1:p
    
    LS(k)= R(k,1,ind);%%%%%%%%%%%%%%R(RULE,DIM,)
    US(k)= R(k,2,ind);
end
CHO_C=[];
CHO_V=[];

for k=1:p
         [c,v,r] = newCENTER2B (numS,US(k),LS(k),k,C(ind,k),p,z);   %%%%%%%%%RD(DIM,SPILT,RANGE)

    
    CHO_C(k,:)=c;
    CHO_V(k,:)=v;
    RD(:,:,k) =r ;
end


for l=1:p

   
    if spilt(l)~=1
    
        CD{1,l}=zeros(spilt(l),p);
        VD{1,l}=zeros(spilt(l),p);
        RRD{1,l}=zeros(p,2,spilt(l));
        for hh=1:spilt(l)
            CD{1,l}(hh,:)=C(ind,:);
            VD{1,l}(hh,:)=V(ind,:);
           
        end
%         VD{1,l}(:,l)=2*V(ind,l)/(spilt(l)+1)
%         VD{1,l}(:,l)=CHO_V(:,l)
        
        if spilt(l)/2-fix(spilt(l)/2)==0
           for hh=1:2:spilt(l)
             
%              CD{1,l}(hh,l)=C(ind,l)-(2*V(ind,l)/(spilt(l)+1) )/2+(hh-1)*2*V(ind,l)/(spilt(l)+1)
%              CD{1,l}(hh+1,l)=C(ind,l)+(2*V(ind,l)/(spilt(l)+1) )/2+(hh-1)*2*V(ind,l)/(spilt(l)+1)
%              
%              CD{1,l}(hh,l)=C(ind,l)-(2*V(ind,l)/(spilt(l)+1) )+(hh-1)*2*V(ind,l)/(spilt(l)+1)
%              CD{1,l}(hh+1,l)=C(ind,l)+(2*V(ind,l)/(spilt(l)+1) )+(hh-1)*2*V(ind,l)/(spilt(l)+1)
%              
             CD{1,l}(hh,l)=CHO_C(l,1);
             CD{1,l}(hh+1,l)=CHO_C(l,2);
           
             VD{1,l}(hh,l)=CHO_V(l,1);
             VD{1,l}(hh+1,l)=CHO_V(l,2);
              
             RRD{1,l}(l,:,hh)=RD(hh,:,l);
             RRD{1,l}(l,:,hh+1)=RD(hh+1,:,l);
             
             
             
           end
           for ll=l+1:p
           
               RRD{1,l}(ll,:,hh)=R(ll,:,ind);
               RRD{1,l}(ll,:,hh+1)=R(ll,:,ind);
           end
           for ll=1:l-1
              RRD{1,l}(ll,:,hh)=R(ll,:,ind);
              RRD{1,l}(ll,:,hh+1)=R(ll,:,ind);
           end 
           
        else
            for hh=1:3:spilt(l)
                
                
                for hhh=0:2
             CD{1,l}(hh+hhh,l)=CHO_C(l,hhh+1);
%              CD{1,l}(hh+1,l)=CHO_C(l,2);
%              CD{1,l}(hh+2,l)=CHO_C(l,3);
%            
             VD{1,l}(hh+hhh,l)=CHO_V(l,hhh+1);
%              VD{1,l}(hh+1,l)=CHO_V(l,2);
%              VD{1,l}(hh+2,l)=CHO_V(l,3);
%              
             RRD{1,l}(l,:,hh+hhh)=RD(hh+hhh,:,l);
%              RRD{1,l}(l,:,hh+1)=RD(hh+1,:,l);
%              RRD{1,l}(l,:,hh+2)=RD(hh+2,:,l);
%             
            end
            end
           for ll=l+1:p
           
                               for hhh=0:2
                                   
               RRD{1,l}(ll,:,hh+hhh)=R(ll,:,ind);
                               end 
           end
           for ll=1:l-1
               for hhh=0:2
              RRD{1,l}(ll,:,hh+hhh)=R(ll,:,ind);
               end 
           end 
           
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if spilt(l)/2-fix(spilt(l)/2)==0
    
    else
       CD{1,l}=C(ind,:);
       VD{1,l}=V(ind,:);
    
    
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if spilt(l)~=1

  end
     
% % % %**************************************************************************
% % % % command between 171   to 233    for each of this  partition  compute loss
% % % % function 
% % % %**************************************************************************
% %  
CCD=cell(1,p);
VVD=cell(1,p);
RRDD=[];
if ind-1==0
   C=[C(ind+1:M,:)];
   V=[V(ind+1:M,:)];
   if M==1&&ind==1
          R=[]
   else
          R=R(:,:,ind+1:M);
   end    
elseif ind==M
    C=[C(1:ind-1,:)];
    V=[V(1:ind-1,:)];
    R=R(:,:,1:ind-1);

else
    C=[C(1:ind-1,:);C(ind+1:M,:)];
    V=[V(1:ind-1,:);V(ind+1:M,:)];
    R1=[];
    R2=[];
    R1(:,:,1:ind-1)=R(:,:,1:ind-1);
    R1(:,:,ind:M-1)=R(:,:,ind+1:M);
    R=R1;
    
end    


for tt=1:p
    CCD{1,tt}=[C;CD{:,tt}];
    VVD{1,tt}=[V;VD{:,tt}];
    
    s=size(CCD{1,tt},1);
    RRDD{1,tt}=[];
    s1=size(C,1);
    if M~=1
    RRDD{1,tt}(:,:,1:s1)=R;
    for ss=s1+1:s
       
        RRDD{1,tt}(:,:,ss)=RRD{1,tt}(:,:,ss-s1);
        

    end    
    else
           
        RRDD{1,tt}=RRD{1,tt};
        

    end    
      
end    
     

for tt=1:p 
    CC=CCD{1,tt};
    VV=VVD{1,tt};
    [LOSS,W]= calLOSSG6(CC,VV,M+spilt(tt)-1);
    lossDUMMY=LOSS;
    WhatDUMMy{1,tt}=W;
    TotalLOSS(tt)=sum(lossDUMMY);
    
end
 
 
 
% % %  %**************************************************************************
% % % % command between 238  compute the best partition that result in minimum loss function  
% % % %**************************************************************************
% [contd ,indd]=min(ll(:,1))
[contd ,indd]=max(TotalLOSS);

% % % STEP  4
% %  %**************************************************************************
% % % command between 243   to  319  add new rule to set of rules
% % %**************************************************************************
C=[CCD{1,indd}];
V=[VVD{1,indd}];
R=[RRDD{1,indd}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      TRAINING THE CENTE AND SIGMA
%%%%%%%%%%%%%%%%%%%%%%% RANDOM SEARCH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
global Law;
Law=size(CCD{1,indd},1);
[LOSS,W]= calLOSSG6(C,V,Law);
WFINAL=W;

% SIG=(1/sqrt(2))./1.5*V
SIG=1./(sqrt(2)*V);
TEMP=ones(N,p);
MEIYU=ones(Law,N);
PHIE=ones(Law,N);
X=[];
for l=1:Law
        
    for k=1:p
            
        TEMP(:,k)= (z(:,k)-C(l,k))*SIG(l,k);

    end    
    TEMP=TEMP.*TEMP;  
    MEIYU(l,:)=sum(TEMP,2)'; 

end 
Z=sum(MEIYU,1);

for l=1:Law
     PHIE(l,:)= MEIYU(l,:) ./ Z;
end

for kk=1:Law
     
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
What=[];
What=inv( X'*X+ alpha* eye(Law*(p+1),Law*(p+1)))*X'*yn;
yhat =X * What;
E=yn - yhat;
Yhat=yhat*( max(y)-min(y))+ min(y);

EE=y-Yhat;
%%%%
% [Imin,Imax,PMAX,PMIN]=PEAKMINMAX(N,yn)
% for i=1:size(Imax)
% EMAX(i)=yn(i)-yhat(i)
% end
% 
% for i=1:size(Imin)
% EMIN(i)=yn(i)-yhat(i)
% end
% %%%%%%%
% AVE(Itteration)=sum((1/N)*(abs(EE)./abs(y)))*100

% perf(Itteration)=mse(E);
% EPEAK(Itteration)=mse(EMIN)*mse(EMAX)*mse(E)
% EPEAK(Itteration)=mse(EMAX)
% EPEAK(Itteration)=sum((1/size(PMAX,2))*(EMAX./PMAX) )*100
% if Itteration>1
%    if EPEAK(Itteration-1)< EPEAK(Itteration)
%       break
%    end    
% 
% end    
% if Itteration>1
%    if AVE(Itteration-1)< AVE(Itteration)
%       break
%    end    
% 
% end    
% AVE(Itteration)=sum((1/N)*(abs(EE)./abs(y)))*100;

perf(Itteration)=mse(E);
% EPEAK(Itteration)=mse(EMIN)*mse(EMAX)*mse(E)
% EPEAK(Itteration)=mse(EMAX)
% EPEAK(Itteration)=sum((1/size(PMAX,2))*(EMAX./PMAX) )*100
% if Itteration>1
%    if EPEAK(Itteration-1)< EPEAK(Itteration)
%       break
%    end    
% 
% end
% [AVEERROR]=CHKVALID(What,Law,C,V)
% AVE(Itteration)=AVEERROR+AVEE(Itteration)

% if Itteration>1
%     if perf(Itteration-1)<perf(Itteration)
%       if AVE(Itteration-1)< AVE(Itteration)
%        break
%       end    
%     end
% end    


% if Itteration>1
%     
%     if AVE(Itteration-1)< AVE(Itteration);
%        break
%       end    
% end    

Itteration=Itteration+1;
end

% Yhat = yhat*((maxU(p+1)-minU(p+1)))+minU(p+1);

UCn=[];
YCn=[];
zc=[];
xc=[];
XC=[];

n = 2;
%%p show the number of input data 
% UT = zeros(NT, n+1);
% UT=[TestP];
% minUT=min(UT);
% maxUT=max(UT);
% for i=1:size(UT,2)
%     PnT(:,i)=(UT(:,i)- minUT(i))./(maxUT(i)-minUT(i));
% end  
UCn=[];
yCn=[];
zc=[];
xc=[];

aanT=aan(N+1:Nt,:)
OUTT=DATA(N+1:Nt,p+1);

UCn=aanT(:,1:p);
% yCn=aanT(:,p+1);
zc=UCn(:,1:p);
xc=UCn(:,1:p);
XC=[];
MEIYUC=[];
PHIEC=[];
yC=AA(:,p+1);
xxc=[];

%use fuzzy system to test data
%use fuzzy system to test data
SigmaT=1./(sqrt(2)*V);
TEMPC=[];
M=Law(1);
for l=1:Law(1)

        
    for k=1:p
            
        TEMPC(:,k)= (zc(:,k)-C(l,k))*SigmaT(l,k);

    end

    TEMPC=TEMPC.*TEMPC;
    MEIYUC(l,:)=sum(TEMPC,2)'; 

end 
ZC=[];
PHIEC=[];
ZC=sum(MEIYUC,1);

for i=1:M
         
        PHIEC(i,:)= MEIYUC(i,:) ./ ZC;
      
end   
XC=[];
%**********************************************************************
% cammand between 91 to 115 for m rule X matrix and What  and yhat matrix compute 
%**********************************************************************
   for i=1:M
     
    
         
        for k=1:p+1
          
          if k==1  
          
              xxc(:,k)=PHIEC(i,:);
        
          else
          
              xxc(:,k) = xc(:,k-1).*PHIEC(i,:)';
          
          end
       
        end  
    
   
XC=[XC,xxc];    
end   
WhatC=WFINAL;
yhatC =XC * WhatC;
toc;
% et=yCn-yhatC;
% perftest=mse(et);


YhatC=yhatC*( max(y)-min(y))+min(y);
ET=YhatC-yC;
PER=mse(ET);
NMSE=(NT*PER)/(mse(yC-mean(yC))*NT);
figure(1);
subplot(1,2,1);
plot(y,'-r');
hold on;
plot(Yhat,'-b');
subplot(1,2,2);
plot(YhatC,'.-b');
hold on ;
plot(yC,'-r');
NRMSE=sqrt(PER)/std(yC);
Corr=corrcoef(YhatC,yC);
