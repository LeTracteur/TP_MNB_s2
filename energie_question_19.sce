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

function energie = cinetique(T,h,n,p)
    //definition des parametres initiaux
    u_0 = zeros(1,n)
    u_m1 = zeros(1,n)
    //nb d'iteration
    N = T/h
    //creation du vecteur qui va contenir les energies cinétiques
    energie = zeros(N/p,63)
    //creation de M
    diago_m = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago_m(i)=m
        end,
    end
    M = diag(diago_m)
    //définition de la diagonale de A et de sa sous diagonale
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
        //on calcul la valeur de l'energie à l'instant khp et on remplit le vecteur
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
//fonction qui trace le graphique de l'énergie cinétique
function plot_cine(energie,h,p)
    xi = [1:63]
    temps = [0:p*h:135]
    clf()
    Sgrayplot(temps,xi,energie)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de l energie cinetique pour h=0.001, p="+string(p), "temps","i")
    xs2jpg(0,"energie_h=0.001,p = "+string(p)+".jpg",1);
endfunction
//on monte la pile au maximum afin d'éviter tout soucis de mémoire
stacksize(268435454)
T = 135
h = 1e-3
n=63
m=1
//On définit p, on prend un diviseur de N = T/h
p = 30

energie=cinetique(T,h,n,p);
plot_cine(energie,h,p);
