
function [linf1,ldiag]=factorise(diago, sous_diag,n)
  //  longueur_diag = size(diago)
  //  n = longueur_diag(2)
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

function y = descente(linf1,ldiag,b,n)
  //  longueur_diag = size(ldiag)
  //  n = longueur_diag(2)
    y = zeros(1,n)
    y(1) = b(1)/ldiag(1)
    for i =2:n
        y(i) = (b(i)- y(i-1)*linf1(i-1))/ldiag(i)
    end
endfunction

function u = remonte(linf1,ldiag,y,n)
 //   longueur_diag = size(ldiag)
 //   n = longueur_diag(2)
    u = zeros(1,n)
    u(n) = y(n)/ldiag(n)
    for i =(n-1):-1:1
        u(i) = (y(i)-linf1(i)*u(i+1))/ldiag(i)
    end
endfunction

function [vitesse] = resolution(T,h,n,m,w)
    //definition des parametres initiaux
    u_0 = zeros(1,n)
    u_m1 = zeros(1,n)
    //nb d'iteration
    N = T/h
    vitesse = zeros(2,N+1);
    //creation de M
    diago_m = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago_m(i)=m
        end,
    end
    M = diag(diago_m)
    //definition des diagonale de A
    sous_diag = (-h^2)*ones(1,n-1)
    
    diago = zeros(1,n)
    diago(1) = M(1,1)+h^2
    for i1 = 2:n
        diago(i1)=M(i1,i1)+2*h^2
    end
    //calcul de la facto de cholesky
    [linf1,ldiag]= factorise(diago, sous_diag)
    for k = 0:N
        t=(k+1)*h
        f1 = tanh(t)*sin(w*t)
        //creation du vecteur b^k
        bk = zeros(1,n)
        bk(1) = 2*M(1,1)*u_0(1) - M(1,1)*u_m1(1) + (h^2)*f1
        for j = 2:n
            bk(j) = 2*M(j,j)*u_0(j) - M(j,j)*u_m1(j)
        end
        vitesse(1,k+1) = (u_0(1)-u_m1(1))/h;
        vitesse(2,k+1) = (u_0(31)-u_m1(31))/h;
        yk = descente(linf1,ldiag,bk);
        ukp1 = remonte(linf1,ldiag,yk);
        u_m1 = u_0;
        u_0 = ukp1;
    end
endfunction

function plot_w()
    T=95
    h=1e-1
    m=0.6
    n=31
    
    ew = zeros(1,13)
    i=0
    for w = [1.2:0.1:2.4]
        i = i+1
        vitesse = resolution(T,h,n,m,w)
        ew(i) = max(abs(vitesse(2,:)))/max(abs(vitesse(1,:)))
    end
    clf()
    w = [1.2:0.1:2.4]
    plot(w,ew,'+')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de e(w) en fonction de w", "w","$e(w)$")
endfunction
plot_w()
