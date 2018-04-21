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

function plot_cine()
    xi = [1:63]
    h=1e-3
    p=65
    T=135
    N = T/h
    temps = [0:65*h:135]
    //definition des parametres initiaux
    u_0 = zeros(1,63)
    u_m1 = zeros(1,63)
    energie = zeros(2077,63)
    M=eye(63,63)
    sous_diag = (-h^2)*ones(1,63-1)
    diago = zeros(1,63)
    diago(1) = M(1,1)+ h^2
    for i1 = 2:63
        diago(i1)=M(i1,i1)+2*h^2
    end
    //calcul de la facto de cholesky
    [linf1,ldiag]= factorise(diago, sous_diag,63)
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
        bk = zeros(1,63)
        bk(1) = 2*M(1,1)*u_0(1) - M(1,1)*u_m1(1) + (h^2)*f1
        for j = 2:63
            bk(j) = 2*M(j,j)*u_0(j) - M(j,j)*u_m1(j)
        end
        if modulo(k,p) == 0 then
            for l = 1:63
                kp = k/65
                energie(kp,l) = 0.5*M(l,l)*((u_0(l)-u_m1(l))/h).^2
            end
        end
        yk = descente(linf1,ldiag,bk,63);
        ukp1 = remonte(linf1,ldiag,yk,63);
        u_m1 = u_0;
        u_0 = ukp1;
    end
    for l = 1:63
        energie(2077,l) = 0.5*M(l,l)*((u_0(l)-u_m1(l))/h).^2
    end
    clf()
    Sgrayplot(temps,xi,energie)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de l energie cinetique pour h=0.001", "temps","i")
    xs2jpg(0,'energie_h=0.001.jpg',1);
endfunction
plot_cine()
