module positions
    use types
    use box_pot_parameters, only: np, sqtwodt, boxl, deltat
    use randomm

    implicit none
    private

    public position, position_ih

contains

    subroutine position_ih(rpos, f, dij, Rz, pbc)
        real(dp), intent(in) :: f(:, :)
        real(dp), intent(inout) :: rpos(:, :)
        real(dp), intent(in) :: dij(:, :)
        real(dp), intent(in) :: Rz(:)
        logical, intent(in):: pbc
        ! Local variables
        integer :: i, j, ij
        integer, parameter :: k = 3 ! Dimensión espacial del sistema
        real(dp) :: temp(np*k)
        real(dp) :: mulout(k)

        temp = 0.0_dp

        do i = 1, k*np, k
            do j = 1, np - 1
                ij = (k*j) - k + 1
                call dgemv('n', k, k, 1.0_dp, dij(i:i + 2, ij:ij + 2), &
                        k, f(:, j), 1, 0.0_dp, mulout, 1)
                temp(i:i + 2) = temp(i:i + 2) + mulout
            end do
        end do

        do i = 1, np
            ij = (k*i) - k + 1
            rpos(1, i) = rpos(1, i) + temp(ij)*deltat + sqtwodt*Rz(ij)
            rpos(2, i) = rpos(2, i) + temp(ij + 1)*deltat + sqtwodt*Rz(ij + 1)
            rpos(3, i) = rpos(3, i) + temp(ij + 2)*deltat + sqtwodt*Rz(ij + 2)
            if (pbc) then
                rpos(:, i) = rpos(:, i) - boxl*dnint(rpos(:, i)/boxl)
            end if
        end do
    end subroutine position_ih

    subroutine position(r, f, pbc)
        real(dp), intent(in) :: f(:, :)
        real(dp), intent(inout) :: r(:, :)
        logical, intent(in):: pbc
        ! Local variables
        integer :: i, j
        integer, parameter :: n = 3 ! Dimensiones espaciales
        real(dp) :: std_rand(n, np)

        ! Distribución normal con desviación estándar igual a sqrt(2*\Delta t)
        do i = 1, n
            do j = 1, np
                std_rand(i, j) = gasdev() * sqtwodt
            end do
        end do

        do i = 1, np
            r(:, i) = r(:, i) + ( f(:, i)*deltat ) + std_rand(:, i)
            if (pbc) then
                r(:, i) = r(:, i) - boxl*dnint(r(:, i)/boxl)
            end if
        end do

    end subroutine position
end module positions
