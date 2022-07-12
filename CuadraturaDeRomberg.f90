program Romberg
use funcion
implicit none 
real*8::a,b,h,suma_i
real*8, allocatable::R(:,:)
integer::i,j,i1,n
write(*,*)'Integracion de numerica con el metodo de Romberg'
write(*,*)'ingrese el intervalo de menor a mayor'
    read(*,*)a,b
write(*,*)'Ingrese J tal que 2^J sean los intevalos equiespaciados'
read(*,*)n
allocate(R(n,n))
h=(b-a)/(2**i)
R(0,0)=(b-a)*0.5*(f(a)+f(b))
do i=1,n 
    suma_i=0
    do i1=1,(2**i-1),2
    suma_i=suma_i+f(a+((b-a)/2**i)*i1)
    end do
    R(i,0)=0.5*(R(i-1,0)+(((b-a)/(2**(i-1)))*suma_i))
    write(*,*)'i=',i,'i0=',i1,'R(i,0)=',R(i,0)
end do
do i=1,n
    do j=i,n
        R(i,j)=((4**j)*R(i,j-1)-R(i-1,j-1))/(4**j-1)
    end do
end do
write(*,*)'La integral mediante el método de Romberg es',R(n,n)
end program Romberg
