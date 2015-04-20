function [f,axis,xticks,yticks] = CTPConstraint(str)
if strcmp(str,'CTP1')
    c1=@(x,y)y-0.858*exp(-0.541*x);
    c2=@(x,y)y-0.728*exp(-0.295*x);
    f={c1,c2};
    axis=[0 1 0.25 1.25];
    xticks=0:0.2:1;
    yticks=0.25:0.25:1.25;
    return
elseif strcmp(str,'CTP2')
    theta = -0.2*pi;
    a = 0.2;
    b=10;
    c=1;
    d=6;
    e=1;
    axis=[0 1 0 2];
    xticks=0:0.1:1;
    yticks=0:0.2:2;
elseif strcmp(str,'CTP3')
    theta = -0.2*pi;
    a = 0.1;
    b=10;
    c=1;
    d=0.5;
    e=1;
    axis=[0 1 0 2];
    xticks=0:0.1:1;
    yticks=0:0.2:2;
elseif strcmp(str,'CTP4')
    theta = -0.2*pi;
    a = 0.75;
    b=10;
    c=1;
    d=0.5;
    e=1;
    axis=[0 1 0 2];
    xticks=0:0.1:1;
    yticks=0:0.2:2;
elseif strcmp(str,'CTP5')
    theta = -0.2*pi;
    a = 0.1;
    b=10;
    c=2;
    d=0.5;
    e=1;
    axis=[0 1 0 2];
    xticks=0:0.1:1;
    yticks=0:0.2:2;
elseif strcmp(str,'CTP6')
    theta = 0.1*pi;
    a = 40;
    b=0.5;
    c=1;
    d=2;
    e=-2;
    axis=[0 1 0 20];
    xticks=0:0.1:1;
    yticks=0:2:20;
elseif strcmp(str,'CTP7')
    theta = -0.05*pi;
    a = 40;
    b=5;
    c=1;
    d=6;
    e=0;
    axis=[0 1 0 2];
    xticks=0:0.1:1;
    yticks=0:0.2:2;
elseif strcmp(str,'CTP8')
    theta = -0.05*pi;
    a = 40;
    b=2.0;
    c=1;
    d=6;
    e=0;
    c1=@(x,y)cos(theta)*(y-e)-sin(theta)*x-a*power(abs(sin(b*pi*power(sin(theta)*(y-e)+cos(theta)*x,c))),d);
    theta = 0.1*pi;
    a = 40;
    b=0.5;
    c=1;
    d=2;
    e=-2;
    c2=@(x,y)cos(theta)*(y-e)-sin(theta)*x-a*power(abs(sin(b*pi*power(sin(theta)*(y-e)+cos(theta)*x,c))),d);
    f={c1,c2};
    axis=[0 1 0 20];
    xticks=0:0.1:1;
    yticks=0:2:20;
    return;
end
c1=@(x,y)cos(theta)*(y-e)-sin(theta)*x-a*power(abs(sin(b*pi*power(sin(theta)*(y-e)+cos(theta)*x,c))),d);
f={c1};
end



