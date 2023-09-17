program approximation
    ! Declarations
    real :: LEFT, RIGHT, h, sigma, eps = 0.000000001
    integer :: N = 20, n0 = 2, m = 2, i, degree
    real, allocatable, dimension(:) :: X, Y, polinom, gradient, vector, tmp
    logical :: isFirstIteration
    LEFT = 1.0
    RIGHT = 2.0
    h = (RIGHT - LEFT) / real(N)

    ! Initialize arrays
    allocate(X(N+1))
    allocate(Y(N+1))

    do i = 1, N + 1
        X(i) = LEFT + real(i - 1) * h
    end do
    write (*, '(a)') 'Welcome. This program is designed to approximate function by polinom.'
    print *
    100 continue
    ! Mode choice
    write(*, '(a)') 'Please, choose mode:'
    write (*, '(a)') 'mode 1 - test function 3.12 * X(i) ** 2 + 2.34 * X(i) + 1.56'
    write(*, '(a)') 'mode 2 - task function'
    read *, i
    print *
    select case(i)
    case(1)
        call getTestFunctionValues(X, Y, N + 1)
        m = 1
    case(2)
        call getTargetFunctionValues(X, Y, N + 1)
    case default
        write(*, '(a)') 'Wrong mode. try again'
        goto 100
        stop
    end select

    allocate(polinom(n0 + m))
    allocate(gradient(n0 + m))
    allocate(vector(n0 + m))
    allocate(tmp(n0 + m))

    ! Initialize arrays
    do i = 1, size(polinom)
        polinom(i) = 0
        gradient(i) = 0
        vector(i) = 0
        tmp(i) = 0
    end do

    do degree = n0, n0 + m - 1
        i = 0
        isFirstIteration = .TRUE.
        sigma = 1000

        do while (sigma >= eps)
            tmp = polinom
            i = i + 1
            call getVector(isFirstIteration, gradient, vector, degree, N, X, Y, polinom)
            polinom = polinom + getLambda(polinom, vector, X, Y, N, degree) * vector

            if (isFirstIteration) then
                isFirstIteration = .FALSE.
            else
                sigma = getSigma(tmp, polinom, degree)
            end if
            write(*,'(a11, i4, a12, f10.4)') 'iteration', i, '   grad(F)=   ', sqrt(dot_product(gradient, gradient))
        end do
        print *
        ! Display results
        write(*, '(a, a11, a11, a9, a6)') ' i', 'X(i)     ', 'Y(i)      ', 'Pol(i)    ', 'diff '
        do i = 1, N + 1
            write(*, '(i2, 3f10.6, f10.6)') i - 1, X(i), &
                    Y(i), getPolinomValue(X(i), polinom, degree), &
                    abs(Y(i) - getPolinomValue(X(i), polinom, degree))
        end do

        print *
        ! Polinom coefficients
        write(*, '(a, i3, a)') 'Polinom coefficients for degree', degree
        do i = 0, degree
            write(*, '(a1, i2, a1, f10.6)') 'a', i, ' = ', polinom(i + 1)
        end do

        print *
        print *
    end do

    deallocate(X, Y, polinom, gradient, vector, tmp)

contains


    real function getPolinomValue(point, array, degree)
        implicit none
        real, intent(in) :: point
        integer :: i
        integer, intent(in) :: degree
        real, dimension(:), intent(in) :: array

        getPolinomValue = 0
        do i = 0, degree
            getPolinomValue = getPolinomValue * point + array(degree + 1 - i)
        end do
    end function getPolinomValue

    subroutine getGradient(degree, N, X, Y, polinom, gradient)
        implicit none
        integer, intent(in) :: degree, N
        real, dimension(:), intent(in) :: X, Y, polinom
        real, dimension(:), intent(out) :: gradient
        integer :: i, j

        do j = 1, degree + 1
            gradient(j) = 0
            do i = 1, N + 1
                gradient(j) = gradient(j) + 2 * (X(i) ** (j - 1)) * (getPolinomValue(X(i), polinom, degree) - Y(i))
            end do
        end do
    end subroutine getGradient

    real function getAverageSquareError(array, X, Y, N, degree)
        implicit none
        real, dimension(:), intent(in) :: array, X, Y
        integer, intent(in) :: N, degree
        integer :: i

        getAverageSquareError = 0
        do i = 1, N + 1
            getAverageSquareError = getAverageSquareError + (Y(i) - getPolinomValue(X(i), array, degree)) ** 2
        end do
    end function getAverageSquareError

    subroutine getVector(isFirstIteration, gradient, vector, degree, N, X, Y, polinom)
        implicit none
        logical, intent(in) :: isFirstIteration
        real, dimension(:), intent(in) :: X, Y, polinom
        real, dimension(:), intent(inout) :: gradient, vector
        integer, intent(in) :: degree, N
        real :: t

        if (isFirstIteration) then
            call getGradient(degree, N, X, Y, polinom, gradient)
            vector = (-1.0) * gradient
        else
            t = dot_product(gradient, gradient)
            call getGradient(degree, N, X, Y, polinom, gradient)
            vector = (-1.0) * gradient + (dot_product(gradient, gradient) / t) * vector
        end if
    end subroutine getVector

    real function getLambda(polinom, vector, X, Y, N, degree)
        implicit none
        real, dimension(:), intent(in) :: polinom, vector
        real, dimension(:), intent(in) :: X, Y
        integer, intent(in) :: N, degree
        real :: a, b, c

        a = getAverageSquareError(polinom - vector, X, Y, N, degree)
        b = getAverageSquareError(polinom, X, Y, N, degree)
        c = getAverageSquareError(polinom + vector, X, Y, N, degree)
        getLambda = (a - c) / (2 * (a - 2 * b + c))
    end function getLambda

    real function getSigma(tmp, polinom, degree)
        implicit none
        real, dimension(:), intent(in) :: polinom
        real, dimension(:) :: tmp
        integer, intent(in) :: degree
        integer :: i
        real :: t

        tmp = tmp - polinom
        getSigma = 1000
        do i = 1, degree + 1
            if (polinom(i) > 0 .or. polinom(i) < 0) then ! polinom(i) != 0
                t = abs(tmp(i) / polinom(i))
                if (t < getSigma) then
                    getSigma = t
                end if
            end if
        end do
    end function getSigma

    ! Функции getTargetFunctionValues и getTestFunctionValues уже определены правильно
    subroutine getTargetFunctionValues(X, Y, N)
        implicit none
        integer :: N, i
        real, dimension(N) :: X, Y

        do i = 1, N
            Y(i)=3 * cos(X(i)) + X(i) * sin(X(i) * X(i))
        end do
    end subroutine getTargetFunctionValues

    subroutine getTestFunctionValues(X, Y, N)
        implicit none
        integer :: N, i
        real, dimension(N) :: X, Y

        do i = 1, N
            Y(i) = 3.12 * X(i) ** 2 + 2.34 * X(i) + 1.56
        end do
    end subroutine getTestFunctionValues
end program approximation