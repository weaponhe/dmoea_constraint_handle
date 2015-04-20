ins='CF';
for i=1:10
    if i>=8
        p = 1035;
        g = 300;
        objs=3;
    else
        p = 600;
        g = 500;
        objs=2;
    end
    
    switch i
        case 1
            lim=[0 0.01];
        case 2
            lim=[0 0.01];
        case 3
            lim=[0 0.4];
        case 4
            lim=[0 0.15];
        case 5
            lim=[0 0.5];
        case 6
            lim=[0 0.2];
        case 7
            lim=[0 0.5];
        case 8
            lim=[0 0.3];
        case 9
            lim=[0 0.08];
        case 10
            lim=[0 1];
    end
    
    
    instance = sprintf('%s%d',ins,i);
    h=figure;
    %PEN
    filepath=sprintf('../LOG/PEN/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    plot(gen,igd,'y*-');
    hold on;
    %CDP
    filepath=sprintf('../LOG/CDP/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    plot(gen,igd,'go-');
    hold on;
    %ADP
    filepath=sprintf('../LOG/ADP/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    plot(gen,igd,'rD-');
    title(instance);
    legend('CMOEA/D-DE-PEN','CMOEA/D-DE-CDP','CMOEA/D-DE-ADP');
    xlabel('Generation');
    ylabel('Average IGD Value');
    ylim(lim);
    picPath = sprintf('./output/%s_LIM.jpg',instance);
    saveas(gcf,picPath);
    close(h);
end