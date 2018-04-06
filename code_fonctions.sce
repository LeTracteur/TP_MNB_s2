function [F]=rhs(Y)
    F(1)=Y(3);
    F(2)=Y(4);
    F(3)=-max((Y(1)-Y(2)),0).^(3/2);
    F(4)=max((Y(1)-Y(2)),0).^(3/2);
endfunction

function [sol]=eulerexp(Y,T,N)
    h = T/N;
    sol = ones(1:N+1);
    for i = 1:N+1
        sol(i) = Y;
        
        Y= Y + h*rhs(Y);
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
        x1(i) = Y(1,i);
        x2(i) = Y(2,i);
        v1(i) = Y(3,i);
        v2(i) = Y(4,i);
    end
    subplot(221)
    plot(x1,temps)
    subplot(222)
    plot(x2,temps)
    subplot(223)
    plot(v1,temps)
    subplot(224)
    plot(v2,temps)
endfunction
