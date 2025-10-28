module pipNN
implicit none 

public 
integer InitState
integer, parameter::II = 16, JJ = 60,KK = 60
double precision, allocatable:: Norm(:)
double precision Enor,Emin

public::II,JJ,KK, Norm
public:: Enor, EMin

contains 


subroutine NNinit(wts)
integer i,j,inp
double precision temp
doubleprecision,allocatable::wts(:)
double precision, allocatable::tmp(:)
double precision  Et

InitState =1 

allocate(wts(II*JJ+JJ*KK+KK+JJ+KK+1))
allocate(Norm(2*II))
    
open(newunit = inp, file ="fort.33")

do i=1,2*II
read(inp,*) Norm(i)
end do 
read(inp,*) Emin
read(inp,*) Enor

close(inp)

open(newunit = inp,file = "weights.energy")

do i=1,II*JJ+JJ*KK+KK+JJ+KK+1
read(inp,*) temp 
wts(i) = temp 
end do 

close(inp)

return
end subroutine 

subroutine NNdeinit(wts)
double precision, allocatable::wts(:)

deallocate(wts,Norm)
InitState=0

return 
end subroutine 

double precision function NNinterpolate(wts,X)
double precision, intent(in):: X(II) 
double precision wts(II*JJ+JJ*KK+KK+JJ+KK+1)
double precision E,temp,  y1(JJ), y2(KK), w1ij,w2jk, w3k,b1j,b2k,b3,b4
integer i,j,k

do j = 1,JJ 
    temp = 0.0d0 
    do i=1, II 
    w1ij = wts((j-1)*II + i)
    temp = temp + w1ij*x(i)
    end do 
    b1j = wts(II*JJ+j)
    y1(j) =tanh(b1j+temp)
end do 
do k=1,KK 
    temp = 0.0d0 
    do j=1,JJ 
        w2jk = wts(II*JJ+JJ+(k-1)*JJ+j)
    temp = temp + w2jk*y1(j)
    end do 
    b2k = wts(II*JJ+JJ+KK*JJ+k)
    y2(k) = tanh(b2k+temp)
end do 
temp = 0.0d0 
do k=1,KK 
    w3k = wts(II*JJ+JJ+KK*JJ+KK+k)
    temp = temp + w3k*y2(k) 
end do 
b3 = wts(II*JJ+JJ+KK*JJ+KK+KK+1)
E = tanh(b3+temp)
NNinterpolate = E
return 
end function 

subroutine PipNorm(X)
double precision X(II)
integer i 

do i=1,II 

X(i) =(X(i)-Norm(2*i-1))*Norm(2*i)-1.0d0

end do 


return 
end subroutine 


subroutine PIPcoord(X,F,position)
double precision, parameter::C1=3.5d0,ph1=1.44386274d0,ph2=2.50084363d0,ph3=2.8877255d0, Ct = 3.0d0, Cr = 1.1d0, Cz = 3.0d0, Cpip = 0.7d0, Cv = 2.0d0 
double precision, parameter::unit= 1.450550d0
double precision X(8),F(16),Z,Fourier(6) ,xx,yy,zz , tx, ty, xn,yn,xb,yb,xh,yh, pipterm(4), fterm(9)
double precision Rhh, Rh1a, Rh2a, Rh1b, Rh2b, Rcb, Rcn, Rch , xa,xc,ya,yc,za,zc, rdamp
character position 

!double precision, parameter:: xn=0.0d0,xb=1.450550d0,xh=2.90110d0,yn=0.0d0,yb=0.0d0,yh=0.0d0
integer i,j,Nx,Ny, symm
!COORDS IN ANGSTROM IN CARTESIAN

Z = 0.5d0*(x(3)+x(6))
fourier = 0.0d0

!Deeply repulsive z or out of bounds correction
if (x(3) .le. 0.75d0 ) x(3) = 0.75d0
if (x(6) .le. 0.75d0 ) x(6) = 0.75d0
if (x(3) .ge. 6.00d0 ) x(3) = 6.00d0
if (x(6) .ge. 6.00d0 ) x(6) = 6.00d0
!Out of bounds site correction
if( x(7) .ge. 1.0d0) x(7) = 1.0d0 
if( x(8) .ge. 1.0d0) x(8) = 1.0d0 
if( x(7) .le. -1.0d0) x(7) = -1.0d0 
if( x(8) .le. -1.0d0) x(8) = -1.0d0 

!print*,x(1)-xb, x(2)-yb, x(4)-xb, x(5)-yb

do i=1,2 
xx =x(1+3*(i-1))
yy =x(2+3*(i-1))
zz =x(3*i)

call closestpt(xx,yy,0,xn,yn)
call closestpt(xx,yy,1,xb,yb)
call closestpt(xx,yy,2,xh,yh)

!print*, 'Distances to N,B,H: ', sqrt((xx-xn)**2+(yy-yn)**2),sqrt((xx-xb)**2+(yy-yb)**2),sqrt((xx-xh)**2+(yy-yh)**2)

Fourier(1+(i-1)*3) = (C1-(cos(ph1*(xx-xn)+ph2*(yy-yn))+cos(ph1*(xx-xn)-ph2*(yy-yn))+cos(ph3*(xx-xn)))+zz**2.0)!*exp(-1.5d0*z)
Fourier(2+(i-1)*3) = (C1-(cos(ph1*(xx-xb)+ph2*(yy-yb))+cos(ph1*(xx-xb)-ph2*(yy-yb))+cos(ph3*(xx-xb)))+zz**2.0)!*exp(-1.5d0*z)
Fourier(3+(i-1)*3) = (C1-(cos(ph1*(xx-xh)+ph2*(yy-yh))+cos(ph1*(xx-xh)-ph2*(yy-yh))+cos(ph3*(xx-xh)))+zz**2.0)!*exp(-1.5d0*z)

end do 

do i=1,6
if(Fourier(i) .lt. 0.0d0) print*, 'error'
end do 

Za = x(7)
Zc = x(8)


if(position .eq. 'O')then 
    xa = 0.0d0
    xc = unit
    ya = 0.0d0
    yc = 0.0d0
    symm = 0
else if (position .eq. 'P') then 
    xa = -0.5d0*unit
    xc = 1.5d0*unit
    ya = sqrt(3.0d0)/2.0d0*unit
    yc = sqrt(3.0d0)/2.0d0*unit
    symm = 0 
else if (position .eq. 'N') then 
    xa = -0.5d0*unit
    xc = unit
    ya = sqrt(3.0d0)/2.0d0*unit
    yc = 0.0d0 
    symm = 1
else if (position .eq. 'B') then 
    xa = 0.0d0
    xc = 1.5d0*unit
    ya = 0.0d0
    yc = sqrt(3.0d0)/2.0d0*unit
    symm = 1
else 
    stop 'Wrong surface site input!'
end if 

Rhh = sqrt((x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2)
Rh1a = sqrt((x(1)- xa)**2 +(x(2)- ya)**2+(x(3)-za)**2)
Rh1b = sqrt((x(1)- xc)**2 +(x(2)- yc)**2+(x(3)-zc)**2)
Rh2a = sqrt((x(4)- xa)**2 +(x(5)- ya)**2+(x(6)-za)**2)
Rh2b = sqrt((x(4)- xc)**2 +(x(5)- yc)**2+(x(6)-zc)**2)

!Deeply repulsive region
if(Rhh .le. 0.3d0) Rhh = 0.3d0 
if(Rh1a .le. 0.75d0) Rh1a = 0.75d0 
if(Rh1b .le. 0.75d0) Rh1b = 0.75d0
if(Rh2a .le. 0.75d0) Rh2a = 0.75d0 
if(Rh2b .le. 0.75d0) Rh2b = 0.75d0
!Asymptotes
if(Rhh .ge. 4.0d0) Rhh = 4.0d0 
if(Rh1a .ge. 5.0d0) Rh1a = 5.00d0 
if(Rh1b .ge. 5.0d0) Rh1b = 5.00d0
if(Rh2a .ge. 5.0d0) Rh2a = 5.00d0 
if(Rh2b .ge. 5.0d0) Rh2b = 5.00d0

xx =0.5d0*(x(1)+x(4))
yy =0.5d0*(x(2)+x(5))
zz =0.5d0*(x(3)+x(6))

Fterm(1) = exp(-zz/Cv)*(Fourier(1)+Fourier(4))**2.0d0
Fterm(2) = exp(-zz/Cv)*(Fourier(2)+Fourier(5))**2.0d0
Fterm(3) = exp(-zz/Cv)*(Fourier(3)+Fourier(6))**2.0d0
Fterm(4) = exp(-zz/Cv)*(Fourier(1)-Fourier(4))**2.0d0
Fterm(5) = exp(-zz/Cv)*(Fourier(2)-Fourier(5))**2.0d0
Fterm(6) = exp(-zz/Cv)*(Fourier(3)-Fourier(6))**2.0d0

call closestpt(xx,yy,0,xn,yn)
call closestpt(xx,yy,1,xb,yb)
call closestpt(xx,yy,2,xh,yh)

Fterm(7) = (C1-(cos(ph1*(xx-xn)+ph2*(yy-yn))+cos(ph1*(xx-xn)-ph2*(yy-yn))+cos(ph3*(xx-xn))))*exp(-zz/Cv)
Fterm(8) = (C1-(cos(ph1*(xx-xb)+ph2*(yy-yb))+cos(ph1*(xx-xb)-ph2*(yy-yb))+cos(ph3*(xx-xb))))*exp(-zz/Cv)
Fterm(9) = (C1-(cos(ph1*(xx-xh)+ph2*(yy-yh))+cos(ph1*(xx-xh)-ph2*(yy-yh))+cos(ph3*(xx-xh))))*exp(-zz/Cv)

if(symm .eq. 0) then 
    pipterm(1) = (exp(-Cpip*rh1a)+exp(-Cpip*rh2a))
    pipterm(2) = (exp(-Cpip*rh1b)+exp(-Cpip*rh2b))
    pipterm(3) = (exp(-Cpip*rh1a)-exp(-Cpip*rh2a))
    pipterm(4) = (exp(-Cpip*rh1b)-exp(-Cpip*rh2b))
else if(symm .eq. 1) then 
    pipterm(1) = rh1a*rh1b+rh2a*rh2b 
    pipterm(2) = (rh1a*rh1b-rh2a*rh2b)**2.0d0
    pipterm(3) = rh1a*rh1b*rh2a*rh2b
    pipterm(4) = (rh1a*rh2b)**2.0d0 + (rh2a*rh2b)**2.0d0
else 
    stop 'wrong symmetry label'
end if 

rdamp = sqrt((xx-0.5d0*(xa+xc))**2.0d0 +(yy-0.5d0*(ya+yc))**2.0d0)

do i=1,9 
    F(i) = Fterm(i)
end do 

do i=1,4 
    F(9+i) = pipterm(i)**2.0d0

end do 
F(14) = Rhh
F(15) = za
F(16) = zc

return
end subroutine 

subroutine closestpt(x,y,mode,xo,yo)
double precision, parameter::unit= 1.450550d0
integer mode
double precision xx,yy,xo,yo , yt,xt ,x,y
integer Nx,Ny 

if (mode .ne. 2) then 

xx=x-mode*unit
yy = y

else if (mode .eq.2) then 
xx = x-0.5d0*unit 
yy = y-sqrt(3.0d0)/2.0d0*unit 
else 
    stop 'mode error'
endif 

xt = xx

if( xx .lt. 0.0d0 ) then 
    xt = xx + real(int(xx/(3.0d0*unit))+1)*3.0d0*unit
    Nx = int(xx/(3.0d0*unit))-1
else

xt = modulo(xt, 3.0d0*unit)
Nx = int(xx/(3.0d0*unit))
end if 

yt =yy

if( yy .lt. 0.0d0 ) then 
    yt = yy + real(int(yy/(sqrt(3.0d0)*unit))+1)*sqrt(3.0d0)*unit
    Ny = int(yy/(sqrt(3.0d0)*unit))-1
else
yt = modulo(yt, sqrt(3.0d0)*unit)
Ny = int(yy/(sqrt(3.0d0)*unit))
end if 

!print*, 'x: ',xx,xt,nx 
!print*,'y: ',yy,yt,ny

if( yt .lt. unit*sqrt(3.00)/2.0d0) then 
    !print*,'bottom'
    if (yt .lt. unit*sqrt(3.0d0)-xt*sqrt(3.0d0)) then
        xt = 0.0d0 
        yt = 0.0d0
    else if(yt .lt. sqrt(3.0d0)*xt - 2.0d0*unit*sqrt(3.0d0)) then 
        xt = 3.0d0*unit
        yt = 0.0d0
    else 
        xt = 1.5d0*unit 
        yt = sqrt(3.0d0)/2.0d0*unit
    end if 

else 
    !print*,'upper'
    if (yt .gt. xt*sqrt(3.0d0)) then
        xt = 0.0d0 
        yt = sqrt(3.0d0)*unit 
    else if(yt .gt. 3.0d0*sqrt(3.0d0)*unit-sqrt(3.0d0)*xt) then 
        xt = 3.0d0*unit
        yt = sqrt(3.0d0)*unit 
    else 
        xt = 1.5d0*unit 
        yt = sqrt(3.0d0)/2.0d0*unit
    end if 

endif 

if (mode .ne. 2) then 
xo = Nx*3.0d0*unit+xt + mode*unit

yo = Ny*sqrt(3.0d0)*unit+yt 

else if (mode .eq. 2) then 
    xo = Nx*3.0d0*unit+xt + 0.5d0*unit

    yo = Ny*sqrt(3.0d0)*unit+yt + sqrt(3.0d0)/2.0d0*unit

else 
    stop 'mode error'
end if 

!print*,'proposed x:', xo,xt

!print*,'proposed y:', yo,yt

return 
end subroutine 

endmodule 