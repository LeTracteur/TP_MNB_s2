//question 15
function [linf1,ldiag]=factorise(diago, sous_diag,n)
    ldiag = zeros(1,n)
    linf1 = zeros(1,n-1)
    //calcul des élements de la diagonale et de la sous diagonale
    ldiag(1) = sqrt(diago(1))
    linf1(1) = sous_diag(1)/ldiag(1)
    for i = 2:(n-1)
        ldiag(i) = sqrt(diago(i)-linf1(i-1).^2)
        linf1(i) = sous_diag(i)/ldiag(i)
    end
    ldiag(n) = sqrt(diago(n)-linf1(n-1).^2)
endfunction

//question 16-17
function y = descente(linf1,ldiag,b,n)
    y = zeros(1,n)
    y(1) = b(1)/ldiag(1)
    for i =2:n
        y(i) = (b(i)- y(i-1)*linf1(i-1))/ldiag(i)
    end
endfunction

//question 16-17
function u = remonte(linf1,ldiag,y,n)
    u = zeros(1,n)
    u(n) = y(n)/ldiag(n)
    for i =(n-1):-1:1
        u(i) = (y(i)-linf1(i)*u(i+1))/ldiag(i)
    end
endfunction

//question 18
function [vitesse,vitesseT4] = resolution(T,h,n,m)
    //definition des parametres initiaux
    u_0 = zeros(1,n)
    u_m1 = zeros(1,n)
    //nb d'iteration
    N = T/h
    //vecteur contenant les valeurs des vitesses 1, 8, 32, 48, 63
    vitesse = zeros(5,N+1);
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
    diago(1) = M(1,1) + h^2
    for i1 = 2:n
        diago(i1)=M(i1,i1)+2*h^2
    end
    //calcul de la facto de cholesky
    [linf1,ldiag]= factorise(diago, sous_diag)
    for k = 0:N
        //calcul de la veleur de f1(t)
        t = (k+1)*h
        if t >=0 & t<=0.5 then
            f1 = t
        elseif t>0.5 & t<=1 then
            f1 = 1 -t
        else
            f1 =0
        end
        //creation du vecteur b^k
        bk = zeros(1,n)
        bk(1) = 2*M(1,1)*u_0(1) - M(1,1)*u_m1(1) + (h^2)*f1
        for j = 2:n
            bk(j) = 2*M(j,j)*u_0(j) - M(j,j)*u_m1(j)
        end
        if t == T/4 then
            //vecteur des vitesses à l'instant t = T/4
            vitesseT4 = ones(1,n) 
            for u = 1:n
                vitesseT4(u)= (u_0(u)-u_m1(u))/h
            end
        end
        vitesse(1,k+1) = (u_0(1)-u_m1(1))/h;
        vitesse(2,k+1) = (u_0(8)-u_m1(8))/h;
        vitesse(3,k+1) = (u_0(32)-u_m1(32))/h;
        vitesse(4,k+1) = (u_0(48)-u_m1(48))/h;
        vitesse(5,k+1) = (u_0(63)-u_m1(63))/h;
        yk = descente(linf1,ldiag,bk);
        ukp1 = remonte(linf1,ldiag,yk);
        u_m1 = u_0;
        u_0 = ukp1;
    end
endfunction

function energie = cinetique(T,h,n)
    //definition des parametres initiaux
    u_0 = zeros(1,n)
    u_m1 = zeros(1,n)
    //nb d'iteration
    N = T/h
    p=65
    energie = zeros(N/p,63)
    
    diago_m = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago_m(i)=m
        end,
    end
    M = diag(diago_m)
    
    sous_diag = (-h^2)*ones(1,n-1)
    diago = zeros(1,n)
    diago(1) = M(1,1) + h^2
    for i1 = 2:n
        diago(i1)=M(i1,i1)+2*h^2
    end
    //calcul de la facto de cholesky
    [linf1,ldiag]= factorise(diago, sous_diag,n)
    for k = 1:N
        t = k*h
        if t >=0 & t<=0.5 then
            f1 = t
        elseif t>0.5 & t<=1 then
            f1 = 1 -t
        else
            f1 =0
        end
        //creation du vecteur b^k
        bk = zeros(1,n)
        bk(1) = 2*M(1,1)*u_0(1) - M(1,1)*u_m1(1) + (h^2)*f1
        for j = 2:n
            bk(j) = 2*M(j,j)*u_0(j) - M(j,j)*u_m1(j)
        end
        if modulo(k,p) == 0 then
            for l = 1:63
                kp = k/p
                energie(kp,l) = 0.5*M(l,l)*((u_0(l)-u_m1(l))/h).^2
            end
        end
        yk = descente(linf1,ldiag,bk,n);
        ukp1 = remonte(linf1,ldiag,yk,n);
        u_m1 = u_0;
        u_0 = ukp1;
    end
    for l = 1:63
        energie(N/p+1,l) = 0.5*M(l,l)*((u_0(l)-u_m1(l))/h).^2
    end
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

function plot_vit(v,h,m)
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

function plot_cine(energie,h)
    xi = [1:63]
    temps = [0:65*h:135]
    clf()
    Sgrayplot(temps,xi,energie)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de l energie cinetique pour h=0.001, p=60", "temps","i")
    xs2jpg(0,'energie_h=0.001.jpg',1);
endfunction

//exemple pour voir si ça marche bien ou pas
//Ad = [1,2,2,2,2]
//Asd = [-1,-1,-1,-1]
//[Lsd,Ld] = factorise(Ad,Asd)
//b = [1,2,3,4,5]
//y1 = descente(Lsd,Ld,b)
//u1 = remonte(Lsd,Ld,y1)

stacksize(268435454)
T = 135
h = 1e-3
n=63
m=1


//[vit,vitT4] = resolution(T,h,n,m);

//plot_vit(vitT4,h,m);
//plot_vitesses1(vit,h,m);
energie=cinetique(T,h,n)
plot_cine(energie,h);

