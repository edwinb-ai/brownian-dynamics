module energies
    use types
    use box_pot_parameters

    implicit none

    private
    public force

contains

    subroutine force(r, f, ener, virial)
        real(dp), intent(in) :: r(:, :)
        real(dp), intent(inout) :: f(:, :)
        real(dp), intent(out) :: ener, virial

        ! Local variables
        integer, parameter :: n = 3 ! Spatial dimension
        integer :: i, j
        real(dp) :: rij, uij
        real(dp) :: frij(n), rposij(n)

        ! Variable intialization
        ener = 0.0_dp
        f = 0.0_dp
        virial = 0.0_dp
        frij(:) = 0.0_dp

        do i = 1, np - 1
            do j = i + 1, np
                uij = 0._dp

                rposij(:) = r(:, i) - r(:, j)
                rposij(:) = rposij(:) - boxl*nint(rposij(:)/boxl)
                rij = norm2(rposij)

                if (rij < rc) then
                    call hardsphere(rij, uij, rposij, frij)

                    ener = ener + uij
                    f(:, i) = f(:, i) + frij
                    f(:, j) = f(:, j) - frij
                    virial = virial + sum(frij * rposij)
                end if
            end do
        end do
    end subroutine force

    subroutine hardsphere(rij, uij, rposij, frij)
        ! Variables de entrada
        real(dp), intent(in) :: rij, rposij(:)
        real(dp), intent(inout) :: frij(:), uij
        ! Variables locales
        real(dp) ::  fij

        if (rij < bpot) then
            uij = (a2/dT)*((1.0_dp/rij)**dlr - (1.0_dp/rij)**dla)
            uij = uij + 1.0_dp/dT
            fij = dlr*(1.0_dp/rij)**(dlr + 1.0_dp) - dla*(1.0_dp/rij)**(dla + 1.0_dp)
            fij = a2*fij/dT
        else
            uij = 0.0_dp
            fij = 0.0_dp
        end if

        frij(:) = fij*rposij(:)/rij
    end subroutine hardsphere
end module energies
