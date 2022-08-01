function sys3 = transfer_func(num,den,Ti,Ts,Tf)
    sys = tf(num, den, Ts);
    sys3 = impulse(sys,Ti:Ts:Tf);
end