updatetype = 'PEN';
PF_Type='POF';
for p=1:10
    PROBLEMS= ['CF1 '; 'CF2 '; 'CF3 '; 'CF4 '; 'CF5 '; 'CF6 '; 'CF7 '; 'CF8 '; 'CF9 '; 'CF10';];
    DIMX    = [10 10 10 10 10 10 10 10 10 10];
    NOP     = [500 500 500 500 500 500 500 5000 5000 5000];
    PROPERTY= ['b.'; 'b.'; 'b.'; 'b.'; 'b.'; 'b.';'b.';'b.';'b.';'b.';];
    [PF,PS] = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );
    
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
    
    if size(PF,1)== 3
        objs = 3;
        pops = 1035;
        gen = 300;
        plot3(PF(1,:),PF(2,:),PF(3,:),deblank(PROPERTY(p,:)),'MarkerSize',2); hold on;
    else
        objs = 2;
        pops = 600;
        gen = 500;
        plot(PF(1,:),PF(2,:),deblank(PROPERTY(p,:)),'MarkerSize',8); hold on;
    end
    set(gca,'FontSize',8);
    xlabel('f1');ylabel('f2');
    title(sprintf('%s Pareto front',deblank(PROBLEMS(p,:))));
    xlim([0 1.2]); ylim([0 1.2]);
    set(gca,'XTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});
    set(gca,'YTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'YTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});
    if size(PF,1)== 3
        set(gca,'ZTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'ZTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});
        view([-45 10]);
    end
    grid off; box on;
   
    set(0,'units','pixel');
    hold on;
    filepath = sprintf('../PF/%s/%s/MOEAD_CF%d(%d)_%d_%d_R1.dat',updatetype,PF_Type,p,objs,pops,gen);
    if size(PF,1)== 3
        [f1,f2,f3] = textread(filepath,'%f %f %f');
        scatter3(f1,f2,f3,'ro');
    else
        [f1,f2] = textread(filepath,'%f %f');
        scatter(f1,f2,'ro');
    end
    legend('Real PF','Obtained PF');
     f = sprintf('output/%s_%s.jpg',deblank(PROBLEMS(p,:)),updatetype);
    saveas(gca,f);
    close(h);
end