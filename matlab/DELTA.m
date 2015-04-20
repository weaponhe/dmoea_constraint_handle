filepath = '..\PF\SA\POF\';
filelist = dir([filepath,'*.dat']);
delta = cell(length(filelist),3);
for i = 1:length(filelist)
    item = regexp(filelist(i).name,'MOEAD_(.*)\((\d)\)_[0-9]*_[0-9]*_R1.dat','tokens');
    ins = item{1}{1};
    nobj = item{1}{2};
    delta{i,1} = ins;
    if strcmp(nobj,'2')
        [f1,f2] = textread([filepath,filelist(i).name],'%f %f');
        N=length(f1);
        d=zeros(1,N-1);
        sum = 0;
        for j=1:N-1
            d(j)=sqrt(pow2(f1(j+1)-f1(j))+pow2(f2(j+1)-f2(j)));
            sum = sum + d(j);
        end
        mean = sum/(N-1);
        s1 = 0;
        for j=1:N-1
            s1 = s1 + abs(d(j)-mean);
        end
        
        delta{i,2} = s1/(N-1)*mean;
        disp(delta);
    end
end

filepath = '..\PF\CDP\POF\';
filelist = dir([filepath,'*.dat']);
for i = 1:length(filelist)
    item = regexp(filelist(i).name,'MOEAD_(.*)\((\d)\)_[0-9]*_[0-9]*_R1.dat','tokens');
    ins = item{1}{1};
    nobj = item{1}{2};
    if strcmp(delta{i,1},ins);
        if strcmp(nobj,'2')
            [f1,f2] = textread([filepath,filelist(i).name],'%f %f');
            N=length(f1);
            d=zeros(1,N-1);
            sum = 0;
            for j=1:N-1
                d(j)=sqrt(pow2(f1(j+1)-f1(j))+pow2(f2(j+1)-f2(j)));
                sum = sum + d(j);
            end
            mean = sum/(N-1);
            s1 = 0;
            for j=1:N-1
                s1 = s1 + abs(d(j)-mean);
            end
            
            delta{i,3} = s1/(N-1)*mean;
            disp(delta);
        end
    end
end


