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
    instance = sprintf('%s%d',ins,i);
    h=figure;
    %PEN
    filepath=sprintf('../LOG/PEN/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    semilogy(gen,igd,'y*-');
    hold on;
    %CDP
    filepath=sprintf('../LOG/CDP/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    semilogy(gen,igd,'go-');
    hold on;
    %ADP
    filepath=sprintf('../LOG/ADP/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);
    [gen,igd] = textread(filepath,'%d	%f');
    semilogy(gen,igd,'rD-');
    title(instance);
    legend('CMOEA/D-DE-PEN','CMOEA/D-DE-CDP','CMOEA/D-DE-ADP');
    xlabel('Generation');
    ylabel('Average IGD Value');
    picPath = sprintf('./output/%s_LOG.jpg',instance);
    saveas(gcf,picPath);
    close(h);
end