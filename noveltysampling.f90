program noveltysampling
use pipNN
implicit none 
double precision x(8),et, x_format(8)
integer i,j,wcount

double precision, allocatable:: wts(:)
double precision, external:: energyget 
character position



read(5,*) position

wcount = II*JJ+JJ*KK+KK+JJ+KK+1


call NNinit(wts)

do i=1, 151 
do j=1,151 

    x(1) = 0.0d0 
    x(2) = 0.0d0

    x(3) = 5.0d0 
    x(4) = 0.74d0 
    x(5) = 90.0d0
    x(6) = 0.0d0 
    x(7) = -1.0d0 + 2.0d0*i/(151.0d0)
    x(8) = -1.0d0 + 2.0d0*j/(151.0d0)

    call ReactToCart(x,x_format)

    Et= energyget(x_format,wts,wcount,position)

    write(14,*),x(7),x(8),Et 


end do   
end do



end program

double precision function energyget(X,wts,wcount,position)
use pipNN
double precision X(8),E, NNx(II)
integer wcount
double precision wts(wcount)
character position


call PIPcoord(X,NNx,position)
call PipNorm(NNx)

E = NNinterpolate(wts,NNx)

E = (E+1.0d0)/Enor+Emin

energyget = E
 
end function
subroutine ReactToCart(ReactCor,CartCor)
double precision ReactCor(8),CartCor(8)
double precision X1,X2,Y1,Y2,Z1,Z2
double precision m1,m2,Rx,Ry,Ex,Ey
double precision X,Y,Z,R,Th,Ph

    X = ReactCor(1)
    Y = ReactCor(2)
    Z = ReactCor(3)
    R = ReactCor(4)
    Th = ReactCor(5)
    Ph = ReactCor(6)

    m1 = 1.0d0 
    m2= 1.0d0

    Ex = x - 0.5d0*y
    Ey = sqrt(3.0d0)/2.0d0*y 

    x = ex 
    y = ey 

    Z1 = Z+r*m1/(m1+m2)*cosd(th)
    Z2 = Z-r*m2/(m1+m2)*cosd(th)

    X1 = X + R*m1/(m1+m2)*sind(th)*cosd(ph)
    X2 = X - R*m2/(m1+m2)*sind(th)*cosd(ph)
    
    Y1 = Y + R*m1/(m1+m2)*sind(th)*sind(ph)
    Y2 = Y - R*m2/(m1+m2)*sind(th)*sind(ph)
    
    CartCor(1) = X1
    CartCor(2) = Y1
    CartCor(3) = Z1
    CartCor(4) = X2
    CartCor(5) = Y2 
    CartCor(6) = Z2 
    CartCor(7) = ReactCor(7)
    CartCor(8) = ReactCor(8)

return
end subroutine