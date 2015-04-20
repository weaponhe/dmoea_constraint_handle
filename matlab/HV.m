ins='CTP';
for i=1:8
    h=figure;
    instance = sprintf('%s%d',ins,i);
    picPath = sprintf('./output/HV_%s.jpg',instance);
    %PEN
    filepath=sprintf('../LOG/PEN/HV/HV_MOEAD-%s(2)-p200-g200.dat',instance);
    [gen,hv] = textread(filepath,'%d	%f');
    plot(gen,hv,'y*-');
    hold on;
    %CDP
    filepath=sprintf('../LOG/CDP/HV/HV_MOEAD-%s(2)-p200-g200.dat',instance);
    [gen,hv] = textread(filepath,'%d	%f');
    plot(gen,hv,'go-');
    hold on;
    %ADP
    filepath=sprintf('../LOG/ADP/HV/HV_MOEAD-%s(2)-p200-g200.dat',instance);
    [gen,hv] = textread(filepath,'%d	%f');
    plot(gen,hv,'rD-');
    le = legend('CMOEA/D-DE-PEN','CMOEA/D-DE-CDP','CMOEA/D-DE-ADP');
    xlabel('Generation');
    ylabel('Average HV Value');
    set(le,'Location','SouthEast');
    title(sprintf('%s HV',instance));
    saveas(gcf,picPath);
    close(h);
end