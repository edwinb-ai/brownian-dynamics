module utils
    use types
    use box_pot_parameters
    implicit none
    private
    public unit_matrix, check_nan, save_timeseries, &
            iniconfig, show_m, save_msd, snapshots, parse_inputs

contains

    subroutine iniconfig(r, d)
        ! Variables de entrada
        real(dp), intent(out) :: r(:,:)
        real(dp), intent(in) :: d
        ! Variables locales
        integer :: i

        ! Se inicializan las posiciones
        r(:,1) = (d / 2.0_dp) - rc

        do i = 2,np
            r(1,i) = r(1,i-1) + d
            r(2,i) = r(2,i-1)
            r(3,i) = r(3,i-1)
            if ( r(1,i) > rc ) then
                r(1,i) = r(1,1)
                r(2,i) = r(2,i-1) + d
                if ( r(2,i) > rc ) then
                    r(1,i) = r(1,1)
                    r(2,i) = r(2,1)
                    r(3,i) = r(3,i-1) + d
                end if
            end if
            
        end do
    end subroutine iniconfig

    subroutine unit_matrix(mat)
        real(dp), intent(inout) :: mat(:, :)
        ! Variables locales
        integer :: i, n

        n = min(size(mat, 1), size(mat, 2))

        mat = 0.0_dp

        forall ( i = 1:n ) mat(i, i) = 1.0_dp

    end subroutine unit_matrix

    subroutine show_m(A)
        real(dp), intent(in) :: A(:,:)
        integer :: ix, n

        n = size(A, 1)

        do ix = 1, n
            write(*,*) A(ix, :)
        end do
    end subroutine show_m

    subroutine check_unity(A)
        real(dp), intent(inout) :: A(:,:)
        integer :: i,n

        n = size(A,1)

        do i = 1, n
            if ( (A(i,i)-1.0_dp) <= 0.000001 ) then
                print*, 'Not unity!'
                stop
            end if
        end do

        print*, 'This matrix is unity!'
    end subroutine check_unity

    subroutine matrix_file(A)
        real(dp), intent(inout) :: A(:,:)
        integer :: i,j,n

        n = size(A,1)

        open(1, file = 'matrix.dat', status='unknown')  
        do i = 1, n
            do j = 1, n
                write(1,'(20G12.4)',advance='no') A(i, j)
            end do
            write(1,*) ''
        end do  
        close(1)
    end subroutine matrix_file

    subroutine check_nan(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        real(dp), intent(in) :: x(:)
        integer :: i,n

        n = size(x)

        do i = 1,n
            if ( ieee_is_nan(x(i)) ) then
                print*, 'Found NaN. Stopping now.'
                stop
            end if
        end do
    end subroutine check_nan

    subroutine save_timeseries(filename,x,y,z)
        real(dp), intent(in) :: x(:,:), y(:,:), z(:,:)
        character(len=*), intent(in) :: filename
        character(len=1024) :: newname
        character(len=8) :: fmt
        character(len=8) :: x1
        integer :: i, j, n, m, u

        n = size(x,1)
        m = size(x,2)
        fmt = '(I5.1)'

        ! Ciclar sobre todas las partículas
        do i = 1, m
            write(x1,fmt) i
            newname = filename//trim(adjustl(x1))//'.dat'
            open(newunit=u, file=newname, status="new")
            do j = 1, n
                write(u,'(f11.8,A,f12.8,A,f12.8)') x(j, i),',',y(j, i),',',z(j, i)
            end do
            close(u)
        end do

    end subroutine save_timeseries

    subroutine save_msd(t, wt, ft, nprom, filename)
        ! Variables de entrada/salida
        real(dp), intent(in) :: t(:), wt(:), ft(:)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nprom

        ! Variables locales
        integer :: u, i

        open(newunit=u, file=filename, status='unknown')
        do i=1,nprom-1
            write(u,'(f15.11,A,f15.11,A,f14.11)') t(i+1),',',wt(i),',',ft(i)
        end do
        close(u)
    end subroutine save_msd

    subroutine snapshots(r, istep, filename)
        ! Arguments
        real(dp), intent(out) :: r(:,:)
        integer, intent(in) :: istep
        character(len=*), intent(in) :: filename

        ! Local variables
        logical :: exists
        integer :: io, i, arrsize

        arrsize = size(r(1, :))

        inquire(file=filename, exist=exists)
        if (exists) then
            open(newunit=io, file=filename, position="append", &
                & status="old", action="write")
            write(io, fmt='(i4)') arrsize
            write(io, fmt='(i6.4)') istep
        else
            open(newunit=io, file=filename, status="new", action="write")
            write(io, fmt='(i4)') arrsize
            write(io, fmt='(i6.4)') istep
        end if

        do i = 1, arrsize
            write(io, fmt='(i1,A,f11.8,A,f12.8,A,f12.8)') 1, ' ', r(1, i), ' ', &
                & r(2, i), ' ', r(3, i)
        end do

        close(io)
    end subroutine snapshots

    subroutine parse_inputs(filename, limt, limg, initconf, avgfreq)
        ! Arguments
        character(len=*), intent(in) :: filename
        integer, intent(out) :: limt, limg, avgfreq
        logical, intent(inout) :: initconf

        ! Local variables
        integer :: u
    
        open(newunit=u, file=filename, status="old", action="read")
        ! Read the inputs from the file
        read(u, *) phi ! Density
        read(u, *) np ! Numer of particles
        read(u, *) with_ih ! Whether to use HI or not
        read(u, *) limt ! Thermalization steps
        read(u, *) limg ! Averaging steps
        read(u, *) avgfreq ! Average frequency for saving values
        read(u, *) initconf
        close(u)
    end subroutine parse_inputs
end module utils