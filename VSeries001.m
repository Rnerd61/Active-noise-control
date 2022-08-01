T = 10;
N = 30;
fs = 5120;
Ts = 1/fs;

Tpp = 10;
Tsp = 10;
t = T*fs + 1;

mu = 10^-1;
%X = sin(2*pi*500*(0:Ts:T));
X = rand(t,1) - 0.5;
%X = downsample(UniversalMil.VarName2,10);
% X = downsample(exp1.VarName2,10);
% X = X/max(X);


PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);
PP = PP(1:40);
SP = SP(1:40);
PP = PP/max(PP);
SP = SP/max(SP);


Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal
e_vseries1 = zeros(t,1);                     %error

Cw1 =  zeros(1, N);
Cw2 =  zeros(N, N);
Cw3 =  zeros(N, N, N);


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
                Cy = Cy + Cw2(i,j)*X(n-i+1)*X(n-j+1);
        end
    end
    
    for i=1:min(n, N)
        for j=i:min(n, N)
            for k=j:min(n, N)
                Cy = Cy + Cw3(i,j,k)*X(n-i+1)*X(n-j+1)*X(n-k+1);
            end
        end
    end
    
    
    Cw_sum=[Cy; Cw_sum(1: end-1)];
    
    Ys(n) = sum(Cw_sum.*SP);
    e_vseries1(n)=Yd(n)+Ys(n);
    
    for i=1:min(n, N)
        temp = 0;
        for m=1:min(n-i,length(SP))
            temp = temp + Shw(m)*X(n-i-m+2);
        end
        Cw1(i) = Cw1(i) - mu*e_vseries1(n)*temp;
   end
                  
    for i=1:min(n, N)
        for j=i:min(n, N)
            temp =0;
            for m=1:min(n-j, length(SP))
                temp = temp +Shw(m)*X(n-i-m+2)*X(n-j-m+2);
            end
            Cw2(i,j) = Cw2(i,j) - mu*e_vseries1(n)*temp;
        end
    end
    
    for i=1:min(n, N)
        for j=i:min(n, N)
            for k=j:min(n, N)
                temp =0;
                for m=1:min(n-k, length(SP))
                    temp = temp +Shw(m)*X(n-i-m+2)*X(n-j-m+2)*X(n-k-m+2);
                end
                Cw3(i,j,k) = Cw3(i,j,k) - mu*e_vseries1(n)*temp;
            end
        end
    end
end
toc;

figure(1);
plot(e_vseries1);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(Yd) 
hold on 
plot(Yd-e_vseries1, 'r');
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal')
hold off


figure(5);
plot(Yd) 
hold on 
plot(Yd-e_vseries1, 'r')
hold on
plot(e_vseries1);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','errror residual')
hold off



function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den, Ts);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end
