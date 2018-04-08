function [F]=rhs(Y)
    F(1)=Y(3);
    F(2)=Y(4);
    F(3)=-max((Y(1)-Y(2)),0).^(3/2);
    F(4)=max((Y(1)-Y(2)),0).^(3/2);
endfunction

function [sol]=eulerexp(Y,T,N)
    h = T/N;
    sol = ones(N+1,4);
    for i = 1:N+1
        sol(i,:) = Y;
        tmp = rhs(Y)
        Y = Y + h*tmp';
    end
endfunction

function plot_euler_explicite(Y,h)
    N = 4/h
    temps = ones(1:N+1);
    for i = 0:N
        temps(i+1) = h*i;
    end
    x1 =ones(1:N+1);
    x2 =ones(1:N+1);
    v1 =ones(1:N+1);
    v2 =ones(1:N+1);
       
    for i = 1:N+1
        x1(i) = Y(i,1);
        x2(i) = Y(i,2);
        v1(i) = Y(i,3);
        v2(i) = Y(i,4);
    end
    subplot(221)
    plot(temps,x1)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de x1 en fonction du temps", "temps","x1(t)")
    
    subplot(222)
    plot(temps,x2)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de x2 en fonction du temps", "temps", "x2(t)")
    
    subplot(223)
    plot(temps,v1)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de V1 en fonction du temps", "temps","v1(t)")
    
    subplot(224)
    plot(temps,v2)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de V2 en fonction du temps", "temps","v2(t)")
    xs2jpg(0,'foo.jpg',1);
endfunction

Y = [0,0,1,0]
u = rhs(Y)
Test = eulerexp(Y,4,40)
plot_euler_explicite(Test,1e-1)
