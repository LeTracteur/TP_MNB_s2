function [linf1,ldiag]=factorise(diago, sous_diag)
    longueur_diag = size(diago)
    n = longueur_diag(2)
    ldiag = zeros(1,n)
    linf1 = zeros(1,n-1)
    //remplissage du vecteur diagonal de L
    ldiag(1) = sqrt(diago(1))
    linf1(1) = sous_diag(1)/ldiag(1)
    for i = 2:(n-1)
        ldiag(i) = sqrt(diago(i)-linf1(i-1).^2)
        linf1(i) = sous_diag(i)/ldiag(i)
    end
    ldiag(n) = sqrt(diago(n)-linf1(n-1).^2)
endfunction

function y = descente(linf1,ldiag,b)
    longueur_diag = size(ldiag)
    n = longueur_diag(2)
    y = zeros(1,n)
    y(1) = b(1)/ldiag(1)
    for i =2:n
        y(i) = (b(i)- y(i-1)*linf1(i-1))/ldiag(i)
    end
endfunction

function u = remonte(linf1,ldiag,y)
    longueur_diag = size(ldiag)
    n = longueur_diag(2)
    u = zeros(1,n)
    u(n) = y(n)/ldiag(n)
    for i =(n-1):-1:1
        u(i) = (y(i)-linf1(i)*u(i+1))/ldiag(i)
    end
endfunction

function [vitesse,vitesseT4,u_0,u_m1] = resolution(T,h,n,m)
    //definition des parametres initiaux
    u_0 = zeros(1,n)
    u_m1 = zeros(1,n)
    //nb d'iteration
    N = T/h
    //creation de M
    diago_m = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago_m(i)=m
        end,
    end
    M = diag(diago_m)
    
    sous_diag = (-h^2)*ones(1,n-1)
    diago = zeros(1,n)
    for i = 1:n
        diago(i)=M(i,i)+2*h^2
    end
    //calcul de la facto de cholesky
    [linf1,ldiag]= factorise(diago, sous_diag)
    vitesse = ones(5,N+1)
    for i = 1:N
        t = i*h
        if t >=0 & t<=0.5 then
            f1 = t
        elseif t>0.5 & t<=1 then
            f1 = 1 -t
        else
            f1 =0
        end,
        //creation du vecteur b^k
        bk = zeros(1,n)
        bk(1) = 2*M(1)*u_0(1) - M(1)*u_m1(1) + (h^2)*f1
        for i = 2:n
            bk(i) = 2*M(i)*u_0(i) - M(i)*u_m1(i)
        end
        if t == T/4 then
            vitesseT4 = ones(1,n) 
            for i = 1:n
                vitesseT4(i)= (u_0(i)-u_m1(i))/h
            end
        end,
        vitesse(1,i) = (u_0(1)-u_m1(1))/h;
        vitesse(2,i) = (u_0(8)-u_m1(8))/h;
        vitesse(3,i) = (u_0(32)-u_m1(32))/h;
        vitesse(4,i) = (u_0(48)-u_m1(48))/h;
        vitesse(5,i) = (u_0(63)-u_m1(63))/h;
        yk = descente(linf1,ldiag,bk)
        ukp1 = remonte(linf1,ldiag,yk)
        u_m1 = u_0
        u_0 = ukp1
    end
    vitesse(1,N+1) = (u_0(1)-u_m1(1))/h;
    vitesse(2,N+1) = (u_0(8)-u_m1(8))/h;
    vitesse(3,N+1) = (u_0(32)-u_m1(32))/h;
    vitesse(4,N+1) = (u_0(48)-u_m1(48))/h;
    vitesse(5,N+1) = (u_0(63)-u_m1(63))/h;
endfunction


function plot_vitesses1(v,h,m)
    N = 135/h
    temps = ones(1,N+1);
    for i = 0:N
        temps(i+1) = h*i;
    end
    clf()
    plot(temps,v(1,:),temps,v(2,:),temps,v(3,:),temps,v(4,:),temps,v(5,:))
    legend('$v_1$','$v_8$','$v_{32}$','$v_{48}$','$v_{63}$')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de la vitesse pour m = "+string(m)+", h=0.001", "temps","v(t)")
    xs2jpg(0,'modele lineaire vitesse_h=0.001.jpg',1);
endfunction

function plot_vi(v,h,m)
    x=[1:63]
    clf()
    plot(x,v, '+')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Vitesse Vi au temps t=T/4", "temps","v(t)")
    xs2jpg(0,'modele lineaire.jpg',1);
endfunction
//exemple pour voir si ça marche bien ou pas
//Ad = [1,2,2,2,2]
//Asd = [-1,-1,-1,-1]
//[Lsd,Ld] = factorise(Ad,Asd)
//b = [1,2,3,4,5]
//y1 = descente(Lsd,Ld,b)
//u1 = remonte(Lsd,Ld,y1)
stacksize(10^8)
T = 135
h = 1e-3
n=63
m=1

[vit,vitT4,a,b] = resolution(T,h,n,m)
//plot_vitesses1(vit,h,m)
//plot_vi(vitT4,h,m)
