program principal
    use types
    use box_pot_parameters
    use tensor, only: ih
    use utils
    use positions, only: position, position_ih
    use energies
    use observables
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

! Variables locales
    real(dp), allocatable :: rpos(:, :), f(:, :), r(:), g(:), h(:)
    real(dp), allocatable :: t(:), wt(:), ft(:)
! Deben ser `allocatable` dado que son arreglos grandes
    real(dp), allocatable :: cfx(:, :), cfy(:, :), cfz(:, :)

    real(dp) :: dr, d
    integer, parameter :: k = 3

    real(dp), allocatable :: dij(:, :)
    real(dp), allocatable :: Rz(:)

    integer :: i, istep, nprom, ncp, s
    real(dp) :: enerpot, epotn, virial, pressure
    logical :: pbc = .true.
    integer :: u, reps
    logical :: init_config, exists
    integer :: therm_steps, prod_steps

    ! Initialize the PRNG
    call random_seed()
    ! Read the input values from a file
    call parse_inputs("input.in", therm_steps, prod_steps, init_config, ncp)
    ! Compute the necessary parameters for the simulation
    rho = 6.0_dp * phi / pi
    boxl = (np/rho)**(1._dp/3._dp)
    rc = boxl/2.0_dp
    d = (1.0_dp/rho)**(1._dp/3._dp)
    dr = rc / mr

    write(output_unit, fmt='(a,1f12.8)') 'The length of the box is: ', boxl
    write(output_unit, fmt='(a,1f12.8)') 'The mean interparticle distance is: ', d
    write(output_unit, fmt='(a,1f12.8)') 'Cut radius: ', rc

    if (with_ih) then
        write(output_unit, fmt='(a)') 'With hydrodynamic interactions'
    else
        write(output_unit, fmt='(a)') 'Without hydrodynamic interactions'
    end if

! Inicializar la memoria de los arreglos
    allocate (rpos(k, np), f(k, np), source=0.0_dp)
    allocate (r(mr), g(mr), h(mr), source=0.0_dp)
    s = np*k
    allocate (dij(s, s), Rz(s), source=0.0_dp)

    ! Check if there is a previous configuration
    inquire(file='finalconBD.dat', exist=exists)
    ! If so, load it
    if (exists .and. init_config) then
        open(newunit=u, file='finalconBD.dat', status="old", action="read")
        write(unit=output_unit, fmt='(a)') 'Reading positions from file...'
        do i = 1, np
            read(u, *) rpos(:, i)
        end do
        close(u)
    ! If not, create a new configuration
    else
        call iniconfig(rpos, d)
    end if
    
    ! Compute initial forces and interactions
    call force(rpos, f, enerpot, virial)
    if (with_ih) call ih(rpos, np, dij, Rz)

!Energy of the initial configuration
    epotn = enerpot/real(np)
    write(output_unit, fmt='(a,f12.8)') 'E/N = ', epotn

! Comienza la termalización
!Periodic boundary conditions; pbc > 0
    open (newunit=u, file='energy_BD.dat', status='unknown')
    do istep = 1, therm_steps
        call position(rpos, f, pbc)
        call force(rpos, f, enerpot, virial)
        epotn = enerpot/real(np)
        pressure = rho + (virial / (3.0 * boxl**3))
        if (mod(istep, 100000) == 0) then
            write(output_unit, fmt='(i8,2f15.8)') istep, epotn, pressure
        end if
        if (mod(istep, therm_steps) == 0) then
            write(u, fmt='(3f16.8)') real(istep), epotn, pressure
        end if
    end do
    close (u)

    write(output_unit, fmt='(a)') 'The system has thermalized'

! Termina la termalización y se guarda la configuración final
    open (newunit=u, file='finalconBD.dat', status='unknown')
    do i = 1, np
        write(unit=u, fmt='(3f16.8)') rpos(:, i)
    end do
    close (u)

    mt = ncp / prod_steps
    nprom = 0
    pbc = .false.
    reps = 1 ! The frequency for updating the RPY tensor for HI

    allocate (cfx(mt, np), cfy(mt, np), cfz(mt, np), source=0.0_dp)
    allocate (t(mt), wt(mt), ft(mt), source=0.0_dp)

    ! Call the value for the HI
    if (with_ih) call ih(rpos, np, dij, Rz)

    ! Loop for production steps
    do i = 1, prod_steps
        if ( with_ih ) then
            call position_ih(rpos, f, dij, Rz, pbc)
        else
            call position(rpos, f, pbc)
        end if
        ! Compute the energy
        call force(rpos, f, enerpot, virial)
        epotn = enerpot/real(np)
        if ( (with_ih) .and. (mod(i,reps) == 0) ) call ih(rpos, np, dij, Rz)

        if (mod(i, 100000) == 0) then
            write(output_unit, fmt='(1i8,1f10.8)') i, epotn
        end if
        if (mod(i, prod_steps) == 0) then
            nprom = nprom + 1
            t(nprom) = deltat*prod_steps*(nprom - 1)
            cfx(nprom, :) = rpos(1, :)
            cfy(nprom, :) = rpos(2, :)
            cfz(nprom, :) = rpos(3, :)
            call gr(rpos, g, dr, .true.)
            ! Save positions to file
            call snapshots(rpos, i, 'production.xyz')
        end if
    end do

    call difusion(nprom, cfx, cfy, cfz, wt, ft)

    call save_msd(t, wt, ft, nprom, 'wt_fself.dat')
    call normalize(g,h,r,dr,nprom,'gr_ih.dat')
end program principal
