T = 1;
N = 20;
K = 2;
Q = 9;
b = 5;

del_X = 1;
fs = 5120;
Ts = 1/fs;
t = 1+ T*fs;

Tpp = 10;
Tsp = 10;

mu_w = 5*10^-9;
mu_s = 10^-9;

X = rand(t,1)-0.5;

PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);

Xs = zeros(N,1);
Xt = zeros(40,1);


Yd = zeros(t,1);
Ys = zeros(t,1);
Yt = zeros(length(SP),1);
e = zeros(t,1);                     

w =  rand(N,1);
q = [zeros(b, 1); reshape(linspace(-1,1,Q), Q, 1); zeros(b,1)];


S = zeros(t, 1);
C = 0.5*[1 -2 1; -2 2 0; 1 1 0];

Shw = SP;
tic;
for n=1:t
    
    for i=1:min(n, length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    
    for i=1:min(n,N)
        S(n) = S(n) + w(i)*X(n-i+1);
    end
    
    u = mod(S(n), 1);
    is = floor(S(n));
    is = mod(is, Q) + 1 + b;
    
    

    qt = q(is:is+2, 1);
    U = [u^2; u; 1];
    Ud = [u; 1; 0];
    
    Yt = [U'*C*qt; Yt(1:end-1,1)];
    Ys(n) = SP'*Yt;
    
    Ysd = Ud'*C*qt;
    
    e(n) = Yd(n) - Ys(n);
    
    for i=1:min(n, N)
        temp = 0;
        for m=1:min(n-i,length(SP))
            temp = temp + Shw(m)*X(n-i-m+2);
        end
        w(i) = w(i) - mu_w*e(n)*temp*Ysd;
   end
    
%    w = w + mu_w*e(n)*Ysd*Xs;
    q(is:is+2) = q(is:is+2) + mu_s*e(n)*C'*U;
end
toc;

figure(1)
plot(Yd)
hold on
plot(Ys)
hold off;


function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den, Ts);
    
    sys3 = impulse(sys,Ti:Ts:Tf);
    sys3 = sys3(1:40);
    sys3 = sys3/max(sys3);

end
