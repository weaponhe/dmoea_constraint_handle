updatetype = 'ADP';
PF_Type='POF';
ins='CTP';
colors=['g','m'];
%CTP8
for i=1:8
    instance = sprintf('%s%d',ins,i);
    picName = sprintf('./output/%s_%s.jpg',instance,updatetype);
    set(0,'units','centimeters');
    position=[0 0 12 10];
    h=figure;
    set(h,'PaperType','A4');
    set(h,'PaperUnits','centimeters');
    set(h,'paperpositionmode','auto');
    set(h,'PaperPosition',position);
    set(h,'units','centimeters');
    set(h,'position',position);
    hold off;
    [f,axis,xticks,yticks]=CTPConstraint(instance);
    for j=1:length(f)
        h=ezplot(f{j},axis);
        set(h,'Color',colors(j));
        hold on
    end
    set(gca,'xtick',xticks,'ytick',yticks);
    
    hold on
    filepath = sprintf('./data/UCTP/PF_DMOEA_%s_UNCONSTRAINT_R0_G200.dat',instance);
    [f3,f4] = textread(filepath,'%f  %f');
    plot(f3,f4,'r.')
    hold off
    hold on
    
    filepath = sprintf('../PF/%s/%s/MOEAD_%s(2)_200_200_R1.dat',updatetype,PF_Type,instance);
    [f1,f2] = textread(filepath,'%f  %f');
    plot(f1,f2,'o');
    
    if i == 1
        legend('Constraint Boundries','Constraint Boundries','Unconstraint PF','Obtained PFs');
    elseif i==8
        legend('Constraint Boundries','Constraint Boundries','Unconstraint PF','Obtained PFs');
    else
        legend('Constraint Boundries','Unconstraint PF','Obtained PFs');
    end
    title(sprintf('%s Pareto Front',instance));
    xlabel('f1');
    ylabel('f2');
    saveas(gcf,picName);
end
close all

