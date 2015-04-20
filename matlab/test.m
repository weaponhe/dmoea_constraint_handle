max_run = 5;

i=1;
ins='CF';
instance = sprintf('%s%d',ins,i);
p = 600;
g = 500;
objs=2;
filepath=sprintf('../LOG/CDP/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat',instance,objs,p,g);;
m=load(filepath);
gen=(0:20:500)';
igd=zeros(26,1);
for j=1:26
    for k=1:max_run
        igd(j)=igd(j)+m(26*(k-1)+j,2);
    end
end
semilogy(gen,igd,'y*-');