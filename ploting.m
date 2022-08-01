function ploting(res)
    
    e_cont = res(:,2);
    Yd = res(:,1);
    e_cont_sq = res(:,3);
    
    figure(1);
    plot(e_cont,'y');
    ylabel('Amplitude');
    xlabel('Discrete time k');
    legend('Noise residue')

    figure(2);
    plot(Yd,'b') 
    hold on 
    plot(Yd-e_cont, 'r')
    hold on
    plot(e_cont,'y');
    ylabel('Amplitude');
    xlabel('Discrete time k');
    legend('Noise signal', 'Control signal','errror residual')
    hold off
    
    figure(3);
    plot(e_cont_sq,'y');
    ylabel('Amplitude');
    xlabel('Discrete time k');
    legend('squared Noise residue')
end