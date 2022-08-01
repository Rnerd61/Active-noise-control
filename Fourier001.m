T = 1;
p = 3;
fs = 5120;
Ts = 1/fs;
N = 30;
t = T*fs+1;

mu = 10^-3.5;

%x = sin(2*pi*500*(0:Ts:T));
%x = downsample(filling.VarName2,10);
x = rand(t,1);

%x = downsample(UniversalMil.VarName2,10);

Tpp = 10;
Tsp = 10;



yd = zeros(t,1);
yw = zeros(t,1);
ys = zeros(t,1);
e_side = zeros(t,1);
e_side_sq = zeros(t,1);

PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);
PP = PP/max(PP);
SP = SP/max(SP);

PP = PP(1:40);
SP = SP(1:40);

SP = [SP; zeros(N+1-length(SP),1)];

Cw0 = zeros(1,1);
Cw1 = zeros(N,1);
Cw2 = zeros(N,N);
Cw3 = zeros(N,N,N);

for n=1:t
    
    %order 0
    yw(n) = yw(n) + Cw0;
    
    %order 1
    for i=1:min(n,N)
        yw(n) = yw(n) + Cw1(i,1)*sin(pi*x(n-i+1)/2);
    end
    
    %order 2
    for i=1:min(n,N)
        for j=i:min(n,N)
            if i==j
                yw(n) = yw(n) + Cw2(i,j)*cos(pi*x(n-j+1));
            else
                yw(n) = yw(n) + Cw2(i,j)*sin(pi*x(n-i+1)/2)*sin(pi*x(n-j+1)/2);
            end
        end
    end
    
    %order 3
    for i=1:min(n,N)
        for j=i:min(n,N)
            for k=j:min(n,N)
                if(i==j && j==k)
                    yw(n) = yw(n) + Cw3(i,j,k)*sin(3*pi*x(n-k+1)/2);
                elseif(i==j)
                    yw(n) = yw(n) + Cw3(i,j,k)*cos(pi*x(n-i+1))*sin(pi*x(n-k+1)/2);
                elseif (j==k)
                    yw(n) = yw(n) + Cw3(i,j,k)*sin(pi*x(n-i+1)/2)*cos(pi*x(n-k+1));
                else
                    yw(n) = yw(n) + Cw3(i,j,k)*sin(pi*x(n-i+1)/2)*sin(pi*x(n-j+1)/2)*sin(pi*x(n-k+1)/2);
                end
            end
        end
    end
    
    for i=1:min(n,length(PP))
        yd(n) = yd(n) + PP(i)*x(n-i+1);
    end
    
    for i=1:min(n,N)
        ys(n) = ys(n) + SP(i)*yw(n-i+1);
    end
    
    e_side(n) = yd(n) + ys(n);
    
    e_side_sq(n) = e_side(n)^2;
    
    %order 0
    Cw0 = Cw0 - mu*e_side(n)*sum(SP);
    
    %order 1
    for i=1:min(n,N)
        temp = 0;
        for m=1:min(n-i,length(SP))
            temp = temp + SP(m)*sin(pi*x(n-i-m+2)/2);
        end
        Cw1(i) = Cw1(i) - mu*e_side(n)*temp;
    end
    
    %order 2
    for i=1:min(n,N)
        for j=i:min(n,N)
            temp = 0;
            for m=1:min(n-j,length(SP))
                if i==j
                    temp = temp + SP(m)*cos(pi*x(n-j-m+2));
                else
                    temp = temp + SP(m)*sin(pi*x(n-i-m+2)/2)*sin(pi*x(n-j-m+2)/2);
                end
            end
            Cw2(i,j) = Cw2(i,j) - mu*e_side(n)*temp;
        end
    end
    
    %order 3
    for i=1:min(n,N)
        for j=i:min(n,N)
            for k=j:min(n,N)
                temp = 0;
                for m=1:min(n-k,length(SP))
                    if(i==j && j==k)
                        temp = temp + SP(m)*sin(3*pi*x(n-i-m+2)/2);
                    elseif(i==j)
                        temp = temp + SP(m)*cos(pi*x(n-i-m+2))*sin(pi*x(n-k-m+2)/2);
                    elseif (j==k)
                        temp = temp + SP(m)*sin(pi*x(n-i-m+2)/2)*cos(pi*x(n-k-m+2));
                    else
                        temp = temp + SP(m)*sin(pi*x(n-i-m+2)/2)*sin(pi*x(n-j-m+2)/2)*sin(pi*x(n-k-m+2)/2);
                    end
                end
                Cw3(i,j,k) = Cw3(i,j,k) - mu*e_side(n)*temp;
            end
        end
    end
    
    disp(n);
end

figure(1);
plot(e_side);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(yd)
hold on
plot(yd-e_side, 'r');
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal')
hold off

figure(3);
plot(yd)
hold on
plot(yd-e_side, 'r')
hold on
plot(e_side);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','errror residual')
hold off

figure(5)
plot(e_side_sq);


function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

sys = tf(num, den, Ts);

sys3 = impulse(sys,Ti:Ts:Tf);

end
