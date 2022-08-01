function res = Vseries(X, PP, SP,N, mu, t)

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal
e_vseries = zeros(t,1);                     %error
e_cont_sq = zeros(t,1); 

Cw1 =  zeros(1, N);
Cw2 =  zeros(N, N);
Cw3 =  zeros(N, N, N);

Xhx1 =  zeros(1, N);
Xhx2 =  zeros(N, N);    
Xhx3 =  zeros(N, N, N); 

Cw_sum = zeros(length(SP), 1);
Shw = SP;

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
                Cy = Cy + Cw2(i,j)*X(n-j+1)*X(n+i-j);
        end
    end
    
    for i=1:min(n, N)
        for j=i:min(n, N)
            for k=j:min(j, N)
                Cy = Cy + Cw3(i,j,k)*X(n-k+1)*X(n+i-k)*X(n+j-k);
            end
        end
    end
    
    
    Cw_sum=[Cy; Cw_sum(1: end-1)];
    
    Ys(n) = sum(Cw_sum.*SP);
    e_vseries(n)=Yd(n)+Ys(n);
    
    temp = 0;
    for i=1:min(n, N)
        temp = temp + Shw(i)*X(n-i+1);
    end
    Xhx1=[temp Xhx1(1:N-1)];
    
    for i=1:min(n, N)
        temp = 0;
        for j=j:min(n, length(SP))
            temp = temp +Shw(j)*X(n-j+1)*X(n+i-j);
        end
        Xhx2(i,:) = [temp Xhx2(i, 1:end-1)];
    end
    
    for i=1:min(n, N)
        for j=i:min(n, N)
            temp=0;
            for k=j:min(n, length(SP))
                temp = temp +Shw(k)*X(n-k+1)*X(n+i-k)*X(n+j-k);
            end
            Xhx3(i,j,:) = cat(3,temp,Xhx3(i,j,1:end-1));
        end
    end
    
   Cw1 = Cw1 - mu*e_vseries(n)*Xhx1;
   Cw2 = Cw2 - mu*e_vseries(n)*Xhx2;
   Cw3 = Cw3 - mu*e_vseries(n)*Xhx3;
    
    e_cont_sq(n) = e_vseries(n)^2;
end

res = [Yd e_vseries e_cont_sq];

end