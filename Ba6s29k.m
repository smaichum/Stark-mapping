clear all; close all;

%%Define variable
ni=29;
Io=42034.902; %Ionization Energy
m=0; %the starting magnetic quantum number
l=ni-1; %max orbital quantum number
H1=zeros(ni); %Prepare to build an expectation value matrix
for i=1:l;
    ep=3/2*ni*sqrt(ni^2-i^2);
    A=(sqrt(((i)^2-m^2)/((2*(i-1)+3)*(2*(i-1)+1))));
    H1(i,i+1)=A*ep;
    H1(i+1,i)=A*ep;
end
H1=H1/5.14e9;
%disp('H1=');disp(H1);
min=0;
max=30;
k=1/3;
h=((max-min)/k);
eb=linspace(min,max,h);
R=109736.88;
I=zeros(ni);
for i=1:ni;
%         if i == 1 
%        I(i,i)=-(Io-41856.50)/2/R;
%        elseif i == 2 
%        I(i,i)=-(Io-41860.981)/2/R;
%        elseif i == 3 
%        I(i,i)=-(Io-41876.831)/2/R;
%        elseif i == 4 
%        I(i,i)=-(Io-41904.486)/2/R;%41904.460
%        elseif i == 5 
%        I(i,i)=-(Io-41903.89)/2/R;
%        elseif i > 5
%        a=i-1;
%        d=0.0585*5^5/(a^5);
%        I(i,i)=-1/2/((ni-d)^2);
%        end
            if i == 1 
        I(i,i)=-(1/(ni-0.1956))^2/2;%-0.21%qdf from 33s1S0 n*=28.8044* 
                                    %34s1S0 n*=29.8048 %0.1956
        elseif i == 2 
        I(i,i)=-(1/(ni+0.0640))^2/2;%+0.091826%qdf from 32p1P1 n*=28.0787%-0.0640
        %33p1P1 n*29.0640 
        elseif i == 3 
        I(i,i)=-(1/(ni+0.3182))^2/2;%+0.3415%qdf from 32d1D2 n*=29.3182
        elseif i == 4 
        I(i,i)=-(1/(ni+0.0075))^2/2;%+0.00754%qdf from 29f1F3 41904.486 n*=29.0075
        elseif i == 5 
        I(i,i)=-(1/(ni-0.0543))^2/2;%-0.0543qdf from 6s29g n*=28.9415 0.0585 <-- only this not good
        %% can't find qdf is 0.0185 in any paper
%         elseif i == 6 
%         I(i,i)=-(1/(ni-0.0185))^2/2;
        elseif i>5 % > 6
        a=i-1;
        d=0.0543*4^5/(a^5);%0.0185*5^5/(a^5);
        I(i,i)=-1/2/((ni-d)^2);
        end
end
%disp(I);
lam=eig(0*H1+I);
for i=1:h;
    H2=2*R*((i*k*H1)+I);
    lm=eig(H2);
    lam=[lam,lm];
end
lam=lam(:,2:end);
lam=lam(1:end-1,:);
% lam=2*R*lam;
% figure(1);hold off;
% plot(eb,lam);
% title('Stark map of nonhydrogen atom n=29') ;
% xlabel('Electric field(V/cm)');ylabel('Energy(cm^-1)')
% lam4=lam(3,:);
% lam=lam.';
% eb=eb.';
% P=[eb,lam];
%dlmwrite('starkmap29.dat',P,'precision','%.16e','delimiter',' ');
%amu*2Rydbergtobe to cm-1
% ep=3/2*29*sqrt(29^2-1^2)
% A=(sqrt(((1)^2-0^2)/((2*(0)+3)*(2*(0)+1))))
% off=A*ep/5.14e9*2*109736.88
% off*200*0.025
%%
%n=32 mapping
Rba=109736.88;
%%Define variable
I=42034.902;
ni2=32;
m2=0;
l2=ni2-1;
H12=zeros(ni2);
for i=1:l2;
    ep2=3/2*ni2*sqrt(ni2^2-i^2);
    A2=(sqrt(((i)^2-(m2)^2)/((2*(i-1)+3)*(2*(i-1)+1))));
    H12(i,i+1)=A2*ep2;
    H12(i+1,i)=A2*ep2;
end
H12=H12/5.14e9;
%disp('H1=');disp(H12);
I2=zeros(ni2);
for i=1:ni2;
     if i == 1 
     I2(i,i)=-(I-41892.94)/2/R;%*(3e8)*(6.63e-34)*(1e2)*(6.242e18);
     elseif i == 2 
     I2(i,i)=-(I-41895.715)/2/R;
     elseif i == 3 
     I2(i,i)=-(I-41907.235)/2/R;%32d n*=29.3241
     elseif i == 4 
     I2(i,i)=-(I-41927.649)/2/R;
     elseif i == 5 
     I2(i,i)=-(I-41927.393)/2/R;
%      elseif i == 6 
%      I2(i,i)=-(1/(ni2-0.0185))^2/2;
     elseif i > 5
     a=i-1;
     d=0.0512*5^5/(a^5);
     I2(i,i)=-1/2/((ni2-d)^2);
     end
end
%disp(I2);
lam2=eig(0*H12+I2);
for i=1:h;
    H22=2*Rba*((i*k*H12)+I2);
    lm=eig(H22);
    lam2=[lam2,lm];
end
lam2=lam2(:,2:end);
lam3=lam2(3,:);
figure(2);hold off ;
plot(eb,lam2);
title('Stark map of nonhydrogen atom n=32') ;
xlabel('Electric field(V/cm)');ylabel('Energy(cm^-1)')
lam2=lam2.';
eb=eb.';
eb=transpose(eb);
%P32=[eb,lam2];
lam1=[lam;lam3];
figure(1);hold off;
plot(eb,lam1);
title('Stark map of nonhydrogen atom n=29') ;
xlabel('Electric field(V/cm)');ylabel('Energy(cm^-1)')
%%
figure(3)
plot(eb,lam,eb,lam2)
title('Stark map of nonhydrogen atom n=29 and n=32') ;
xlabel('Electric field(V/cm)');ylabel('Energy(cm^-1)')
figure(4)
lam4=[lam;lam3];
plot(eb,lam4)
title('Simulated Energy level of Rydberg Barium atom in Stark effect around n=29 ') ;
xlabel('Electric field(V/cm)');ylabel('Energy(cm-1)')
nn=size(lam4);
nn=nn(1);
llam=[];
for i=1:nn
    llam=[llam;lam3];
end
lamm=lam4-llam;
%lamm=lamm*1e2;%m-1
% dP=dP/8065.54;%eV
% dP=dP/6.6e-7/2/pi;%GHz
%lamm=3e8*lamm/1e9;
lamm=lamm*30;
kk=lamm(1,end);
lamm=lamm-kk+45;
figure(5)
plot(eb,lamm)
title('Simulated Transition Energy of Rydberg Barium atom in Stark effect around n=29 ') ;
xlabel('Electric field(V/cm)');ylabel('Frequency(GHz)')
A = importdata('stark_Aj.txt');
A=A.data;
f=A(:,1)-kk+45;
e=A(:,2:end);
% figure(7)
% plot(e+25,flip(f))
% title('Observation Result ') ;
% xlabel('Electric field(V/cm)');ylabel('Frequency(GHz)')
e=e+20.+0.6;
ff=flip(f)-200+40-9-0.25-22-1.7;
% %9th line
f9=flip(f)-200+40-9-0.25-22-0.;
e9=e(:,9)-5.8;
% %8th line
f8=flip(f)-200+40-9-0.25-22-1.1;
e8=e(:,8)-4.;
% %7th line
f7=flip(f)-200+40-9-0.25-22-1.4;
e7=e(:,7)-2.6;
% %6th line
f6=flip(f)-200+40-9-0.25-22-1.6;
e6=e(:,6)-1.4;
% %5th line
f5=flip(f)-200+40-9-0.25-22-1.7;
e5=e(:,5);
% %4th line
e4=e(:,4)+1.17;
f4=flip(f)-200+40-9-0.25-22-1.6;
%3rd line
e3=e(:,3)+2.9+0.1;%e(:,3)+2.9+0.
f3=flip(f)-200+40-9-0.25-22-2.1;
%2nd line
e2=e(:,2)+2.9+2.4-0.3;
f2=flip(f)-200+40-9-0.25-22-2.4;
%1st line
e1=e(:,1)+2.9+2.4+1.2;
f1=flip(f)-200+40-9-0.25-22-3.4;
%f=flip(f)-199.0;
%f=(flip(f)-235.0);
figure(8)
plot(eb,lamm,e5,f5,e4,f4,e3,f3,e2,f2,e1,f1,e6,f6,e7,f7,e8,f8,e9,f9)
axis([0 20 70 125])
title('Transition Energy of Rydberg Barium atom in Stark effect around n=29 with obervation') ;
xlabel('Electric field(V/cm)');ylabel('Frequency(GHz)')
figure(9)
plot(eb,lamm,e5,f5,e4,f4,e3,f3,e2,f2,e1,f1,e6,f6,e7,f7,e8,f8,e9,f9)
title('Inside Transition Energy of Rydberg Barium atom in Stark effect around n=29 with obervation') ;
xlabel('Electric field(V/cm)');ylabel('Frequency(GHz)')
% axis([0 15 -110 -70])
figure(7)
plot(e5,f5,e4,f4,e3,f3,e2,f2,e1,f1,e6,f6,e7,f7,e8,f8,e9,f9)
title('Observation Result ') ;
xlabel('Electric field(V/cm)');ylabel('Frequency(GHz)')
figure(6)
plot(eb,lamm,e,ff)
set(gcf,'position',[00 00 500 500])
axis([-10 30 70 125])
%spacing is around 1.5
% AA=A-[A(2:end,:);0.*A(1,:)];
lamA=[lamm(:,2:end),0.*lamm(:,end)];
lammsp=lamm'-lamA';
% E=[e1,e2,e3,e4,e5,e6,e7,e8,e9,f1,f2,f3,f4,f5,f6,f7,f8,f9];
% dlmwrite('starkmap29to32d_observe.txt',E,'precision','%.16e','delimiter',' ');
% dlmwrite('starkmap29to32d_simulation.txt',[eb,lamm],'precision','%.16e','delimiter',' ');