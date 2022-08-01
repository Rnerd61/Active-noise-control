function res = chebyshev(x, PP, SP, N1,p, mu, t)

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal

Cw=zeros(N1, p);                          %weights
e_cont=zeros(t,1);                    %error
e_cont_sq =zeros(t,1); 


Xhx=zeros(N1, p);

Xw =  zeros(t, 1);

for n=1:t
    
    for i=1:min(n,length(PP))
        Yd(n) = Yd(n) + PP(i)*x(n-i +1);
    end
    
    for i=1:p
        for j=1:min(n,N1)
            Xw(n) = Xw(n) + Cw(j,i)*T_func(x(n-j+1), i);
        end
    end
    
    
    for i=1:min(n,N1+1)
        Ys(n) = Ys(n) + SP(i)*Xw(n-i+1);
    end
    
    e_cont(n)=Yd(n)+Ys(n);

    for i=1:p
        temp =0;
        for j=1:min(n,length(SP))
            temp = temp + SP(j)*T_func(x(n-j+1), i);
        end
        Xhx(:,i)= [temp; Xhx(1:end-1,i)];
    end

    for i=1:p
        Cw(:,i) = Cw(:,i) - mu*e_cont(n)*Xhx(:,i);
    end
    
    e_cont_sq(n) = e_cont(n)^2;
end

res = [Yd e_cont e_cont_sq];
end

function res = T_func(x, n)
    if n<=1
        res = x^n;
    else
        res = 2*x*T_func(x, n-1) - T_func(x, n-2);
    end
end