function res = Cox_deBoor(i, p,u,Q)
if p==0
    if u>=Q(i) && u<Q(i)
        res =1;
    else
        res =0;
    end
elseif p==1
    t1 = (u - Q(i))/(Q(i+p) - Q(i));
    t2 = (Q(i+p+1) - u)/(Q(i+p+1) - Q(i+1));
    if u>=Q(i) && u<Q(i+1)
        res= t1;
    elseif u>=Q(i+1) && u<Q(i+2)
        res = t2;
    else
        res =  0;
    end
else
    t1 = (u - Q(i))/(Q(i+p) - Q(i));
    t2 = (Q(i+p+1) - u)/(Q(i+p+1) - Q(i+1));
    
    res = t1*Cox_deBoor(i, p-1,u,Q) + t2*Cox_deBoor(i+1, p-1,u,Q);
end

end