function res = Hamerstein(X, PP, SP, N1,p, mu,t)

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t, 1);                      %Control Signal
e_cont2 = zeros(t,1);                    %error
e_cont_sq = zeros(t,1);

Cw=zeros(N1, p);                          %weights
Xw = zeros(N1, p);
 
Cw_sum = zeros(length(SP), 1);

for n=1:t
    
    for i=1:min(n,length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    
    Cy = 0;
    for i=1:p
        for j=1:min(n, N1)
               Cy = Cy + Cw(j,i)*(X(n-j+1)^i);
        end
    end
    
    
    Cw_sum=[Cy; Cw_sum(1: end-1)];             
    Ys(n) = sum(Cw_sum.*SP);
    
    e_cont2(n)=Yd(n)+Ys(n);

    for i=1:p
        temp =0;
        for j=1:min(n,length(SP))
            temp = temp + SP(j)*(X(n-j+1)^i);
        end
        Xw(:,i)= [temp; Xw(1:end-1,i)];
    end
    
    for i=1:p
        Cw(:,i) = Cw(:,i) - mu*e_cont2(n)*Xw(:,i);
    end
    
    e_cont_sq(n) = e_cont2(n)^n;
end

res = [Yd e_cont2 e_cont_sq];
end