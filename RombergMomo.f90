program ques_03
    implicit none
    integer :: i, j, k, n, l
    real :: a, b, exact, x(0:100), y(0:100), h, tr, rom(20,20)

    exact = 0.0887553; n = 5!n = 5 decimal places
    a = 0.0; b = 3.1416/4! a & b are lower and upper limit respectively

    do i = 1, n
        l = 2**(i-1)
        h = (b-a)/l
        do j = 0, l
            x(j) = j*h
            y(j) = (x(j))**2*sin(x(j))
        end do
        !Trapezoidal rule for getting 1st column of Romberg equation
        tr = 0
        do j = 1, l-1
            tr = tr + y(j)
        end do
        rom(i,1) = h*(y(0)+2*tr+y(l))/2.
    end do

    do j = 2, n
        do i = j, n
            rom(i,j)=rom(i,j-1)+(rom(i,j-1)-rom(i-1,j-1))/(4**(j-1)-1) !!Romberg Theorem
        end do
    end do

    write(*,*)"Extrapolation table are shown bellow: "
    write(*,*)
    write(*,*)"   O(h2)   ","    O(h4)   ","    O(h6)   ","    O(h8)   ","    O(h10)   "
    write(*,*)" --------  ","  --------  ","  --------  ","  --------  ","  ---------  "

    do i=1,n
        write(*,'(5f12.8)')(rom(i,j),j=1,i)  !!For table
    end do

    write(*,*)
    write(*,*)"Exact result is                          : ",exact
    write(*,*)"Result by Romberg's Integration method is: ",rom(n,n)
    write(*,*)"Absolute Error is                        : ",abs(exact-rom(n,n))

end program