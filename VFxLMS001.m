T = 10;
N = 30;
fs = 5120;
Ts = 1/fs;
t = 1+ T*fs;

Tpp = 10;
Tsp = 10;

mu = 10^-3;
%X = sin(2*pi*500*(0:Ts:T));
X = rand(t,1) - 0.5;
% X = downsample(exp1.VarName2,10);
% X  = X/max(X);

PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);

PP = PP(1:40);
SP = SP(1:40);

PP = PP/max(PP);
SP = SP/max(SP);

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal
e_vfxlms = zeros(t,1);                     %error

Cw1 =  zeros(1, N);
Cw2 =  zeros(N, N);

Xhx1 =  zeros(1, N);
Xhx2 =  zeros(N, N);

Cw_sum = zeros(length(SP), 1);
Shw = SP;

tic;
for n=1:t
    
    for i=1:min(n, length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    
    Cy = 0;
    for i=1:min(n,N)
        Cy = Cy + Cw1(i)*X(n-i+1);
    end

    for i=1:min(n, N)
        for j=i:min(n, N)
                Cy = Cy + Cw2(i,j)*X(n-j+1)*X(n-i+1);
        end
    end
  
    Cw_sum=[Cy; Cw_sum(1: end-1)];
    
    Ys(n) = sum(Cw_sum.*SP);
    e_vfxlms(n)=Yd(n)+Ys(n);
    
    Xhx1 = [0 Xhx1(1,1:end-1)];
    
    for m=1:min(n,length(SP))
        Xhx1(1) = Xhx1(1) + Shw(m)*X(n-m+1);
    end
              
    
    Xhx2 = [zeros(1,N); Xhx2(1:end-1,:)];
    Xhx2 = [zeros(N,1) Xhx2(:, 1:end-1)];
    
    for j=1:min(n, N)
        for m=1:min([n-j, length(SP)])
            Xhx2(1,j) = Xhx2(1,j) +Shw(m)*X(n-m+1)*X(n-j-m+2);
        end
    end
    
    Cw1 = Cw1 - mu*e_vfxlms(n)*Xhx1;
    Cw2 = Cw2 - mu*e_vfxlms(n)*Xhx2;
    
end
toc;

figure(1);
plot(e_vfxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(Yd) 
hold on 
plot(Yd-e_vfxlms, 'r');
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal')
hold off


figure(5);
plot(Yd) 
hold on 
plot(Yd-e_vfxlms, 'r')
hold on
plot(e_vfxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','errror residual')
hold off



function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den, Ts);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end
