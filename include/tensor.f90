module tensor
use types
use outerprod
use utils, only: unit_matrix
use randomm, only: gasdev
use ieee_arithmetic, only: ieee_is_nan

implicit none
private
public ih

contains

    subroutine ih(rpos, k, mat_a, R)
        ! Varibles de entrada y salida
        integer, intent(in) :: k   ! numero de submatrices (particulas)
        real(dp), intent(in) :: rpos(:,:)
        real(dp), intent(out) :: R(:)
        real(dp), intent(out) :: mat_a(:, :)

        ! Variables locales
        integer, parameter :: n = 3 !la dimension de la matriz (x,y,z)
        integer :: s
        real(dp) :: Dij(n, n)
        real(dp) :: part1(n)
        real(dp) :: part2(n)
        real(dp) :: ident(n, n)
        real(dp), allocatable :: Xr(:)
        integer :: ix, ij, pos, il   ! indices
        real(dp), allocatable :: sigma(:,:)

        ! Tamaño total de elementos, partícula por cada dimensión espacial
        s = n*k
        allocate( sigma(s,s), Xr(s) )

        ! Inicialización de los arreglos
        call unit_matrix( ident )
        R = 0.0_dp
        Xr = 0.0_dp
        mat_a = 0.0_dp
        Dij = 0.0_dp

        ! Distribucion normal estándar
        do ix = 1,s
            Xr(ix) = gasdev()
        end do

        ! Hacer una matriz identidad a bloques
        do ix = 1, s, n
            mat_a(ix:ix+n-1,ix:ix+n-1) = ident
        end do

        do ix = 1, k
            il = (ix * n) - n + 1
            do ij = ix+1, k
                pos = (ij*n) - n + 1

                part1 = rpos(:, ix)
                part2 = rpos(:, ij)

                call matrizDij( part1,part2,Dij,n )

                mat_a(il:il+n-1,pos:pos+n-1) = Dij
                mat_a(pos:pos+n-1,il:il+n-1) = Dij
            end do
        end do

        ! Se añade un pequeño valor para hacer a la matriz
        ! estable numéricamente
        sigma = mat_a + 0.00001_dp
        !! # Se escoge cualquiera de los dos
        ! call cholesky_method( sigma,s,Xr,R )
        call krylov_method( sigma,s,Xr,R )

        ! Liberar memoria
        deallocate( sigma,Xr )
    end subroutine ih

    subroutine matrizDij(part1, part2, Dij, n)
        ! Variables de entrada y salida
        integer, intent(in) :: n
        real(dp), intent(in) :: part1(:)
        real(dp), intent(in) :: part2(:)
        real(dp), intent(out) :: Dij(:,:)

        ! Variables locales
        real(dp) :: rij,sqrddist,xij,yij,zij
        real(dp) :: ident(n, n), prodout(n, n)
        real(dp) :: temp(n)
        real(dp), external :: dnrm2 ! Norma euclidiana de BLAS

        ! Calcular las distancias
        xij = part1(1) - part2(1)
        yij = part1(2) - part2(2)
        zij = part1(3) - part2(3)
        ! rij = xij**2 + yij**2 + zij**2
        temp = [xij, yij, zij]
        sqrddist = dnrm2( 3,temp,1 )
        rij = sqrddist ** 2

        ! calculando Dij
        Dij = 0.0_dp
        prodout = 0.0_dp

        call unit_matrix( ident )

        if (sqrddist >= 1.0_dp) then
            prodout = outerprod_d( temp, temp ) / rij
            Dij = ident + prodout
            
            Dij = Dij + ( (ident/3.0_dp) - prodout )/(2.0_dp*rij)
            
            Dij = 3.0_dp*Dij/(8.0_dp*sqrddist)
        else
            Dij = (1.0_dp-(9.0_dp*sqrddist/16.0_dp))*ident
            Dij = Dij + 3.0_dp*outerprod_d( temp, temp )/(16.0_dp*sqrddist)
        end if
    end subroutine matrizDij

    subroutine cholesky_method(sigma,s,xr,r)
        !  Variables de entrada
        real(dp), intent(in) :: sigma(:,:),xr(:)
        real(dp), intent(inout) :: r(:)
        integer, intent(in) :: s

        ! Variables locales
        integer :: info
        character(1) :: uplo
        real(dp), allocatable :: temp(:,:)
        allocate( temp(s,s) )

        ! Descomposición de Cholesky
        uplo = 'L'
        call dpotrf( uplo,s,sigma,s,info )
        if ( info .ne. 0 ) print*, "No se puede descomponer"
        ! Copiar la descomposición a una matriz nueva
        call dlacpy( uplo, s, s, sigma, s, temp, s )
        ! Multiplicar L*Xr para obtener el vector de números aleatorios
        call dgemv( 'n',s,s,1.0_dp,temp,s,Xr,1,0.0_dp,R,1 )

        deallocate(temp)
    end subroutine cholesky_method

    subroutine krylov_method(sigma,s,xr,r)
        !  Variables de entrada
        real(dp), intent(in) :: sigma(:,:),xr(:)
        real(dp), intent(inout) :: r(:)
        integer, intent(in) :: s
        ! Variables locales
        real(dp), allocatable :: w(:),eid(:),v(:),res(:)
        real(dp), allocatable :: temp(:,:),vm(:,:),h(:,:)
        integer :: n, m, i, k
        real(dp), external :: dnrm2,ddot ! Operaciones de BLAS
        real(dp) :: znorm, dump

        if ( (size(r) .ne. s) .or. (size(Xr) .ne. s) ) stop 'Tamaños incorrectos'

        ! Inicializar arreglos
        m = 50 ! Número de pasos de Lanczos
        allocate( vm(s,m),h(m,m),w(s),eid(m) )
        allocate( temp(m,m),v(m),res(s) )
        vm(:,:) = 0.0_dp
        h(:,:) = 0.0_dp
        w(:) = 0.0_dp
        v(:) = 0.0_dp
        r(:) = 0.0_dp
        res(:) = 0.0_dp
        temp(:,:) = 0.0_dp
        ! eid es el primer vector de la matrix identidad
        eid(:) = 0.0_dp
        eid(1) = 1.0_dp

        ! Calcular el primer vector base
        znorm = dnrm2( s,xr,1 )
        vm(:,1) = xr / znorm
        ! Comienzan las iteraciones de Lanczos
        do i = 1,m
            call dgemv( 'n',s,s,1.0_dp,sigma,s,vm(:,i),1,0.0_dp,w,1 )
            if ( i > 1 ) then
                w = w - ( h(i-1,i) * vm(:,i-1) )
            end if
            k = size(w, 1)
            h(i,i) = ddot( k,w,1,vm(:,i),1 )
            if ( i < m ) then
                w = w - ( h(i,i) * vm(:,i) )
                k = size(w, 1)
                h(i+1,i) = dnrm2( k,w,1 )
                h(i,i+1) = h(i+1,i)
                vm(:,i+1) = w / h(i+1,i)
            end if
        end do
        call sqrt_matrix( h,temp )
        call dgemv( 'N',m,m,1.0_dp,temp,m,eid,1,0.0_dp,v,1 )
        call dgemv( 'N',s,m,1.0_dp,vm,s,v,1,0.0_dp,res,1 )
        r = znorm * res

        ! Se libera la memoria
        deallocate( v,vm,h,w,temp,eid,res )
    end subroutine krylov_method

    subroutine sqrt_matrix(a,b)
        ! Variables de entrada
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: b(:,:)
        ! Variables locales
        integer :: lda,n,i,nb,info,lwork
        real(dp), allocatable :: w(:),work(:)
        real(dp), allocatable :: diagw(:,:),temp(:,:)

        lda = size(a, 1)
        n = size(a, 2)
        nb = 64
        lwork = nb * n
        allocate( w(lda),diagw(lda,n),work(lwork) )
        allocate( temp(lda,n) )
        diagw = 0.0_dp
        b = 0.0_dp
        w = 0.0_dp
        work = 0.0_dp
        
        ! Copiar los datos de `a` en `b`
        call dlacpy( 'All',lda,n,a,lda,b,lda )

        ! Descomposición espectral de `b`
        call dsyev( 'V','U',n,b,lda,w,work,lwork,info )
        if (info .ne. 0) stop 'No se pudo descomponer espectralmente'

        ! Construir la matriz diagonal, habiendo calculado la raíz cuadrada
        w = dsqrt(w)
        forall ( i = 1:lda ) diagw(i,i) = w(i)

        ! Reconstruir la matriz
        ! Primero, Lambda * U**T
        call dgemm( 'N','T',lda,n,n,1.0_dp,diagw,lda,b,lda,0.0_dp,temp,lda )
        ! Luego, U * temp = U * (Lambda * U**T)
        call dgemm ( 'N','N',lda,n,n,1.0_dp,b,lda,temp,lda,0.0_dp,b,lda )

        ! Liberar la memoria
        deallocate( w,diagw,work,temp )
    end subroutine sqrt_matrix
end module tensor