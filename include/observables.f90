module observables
    use types
    use box_pot_parameters
    implicit none
    save
    public gr, normalize, difusion
contains
    subroutine gr(r,g,dr,pbc)
        real(dp), intent(in) :: r(:,:)
        real(dp), intent(inout) :: g(:)
        real(dp), intent(in) :: dr
        logical, intent(in) :: pbc
        real(dp), external :: dnrm2
        
        ! Local variables
        integer :: i, j, nbin
        real(dp) :: rposij(3), rij

        do i = 1,np-1
            do j = i+1,np
                rposij = r(:,i) - r(:,j)
                if (pbc) then
                    rposij = rposij - boxl*dnint(rposij/boxl)
                end if
                rij = dnrm2( 3,rposij,1 )
                if (rij < rc) then
                    nbin = idnint(rij/dr) + 1
                    if (nbin <= mr) then
                        g(nbin) = g(nbin) + 2.0_dp
                    end if
                end if
            end do
        end do
    end subroutine gr

    subroutine normalize(g,h,r,dr,nprom,filename)
        ! Variables de entrada/salida
        real(dp), intent(inout) :: g(:), h(:), r(:)
        real(dp), intent(in) :: dr
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nprom

        ! Variables locales
        integer :: u, i
        real(dp) :: dv, fnorm, graux, hraux

        open(newunit=u, file=filename, status='unknown')
        write(u,'(3f16.8)') r(1), g(1)

        do i=2,mr
            r(i)=(i-1)*dr
            dv=4.0_dp*pi*r(i)*r(i)*dr
            fnorm=boxl**3.0_dp/( np**2.0_dp * nprom*dv )
            graux=g(i)*fnorm
            hraux=graux-1.0_dp
            g(i)=graux
            h(i)=hraux
            write(u,'(3f16.8)')r(i),graux,hraux
        end do
        close(u)
    end subroutine normalize

    ! This sobroutine computes the mean-square displacement
    subroutine difusion(nprom, cfx, cfy, cfz, wt, ft)
        real(dp), intent(in) :: cfx(:,:), cfy(:,:), cfz(:,:)
        real(dp), intent(inout) :: wt(:), ft(:)
        integer, intent(in):: nprom
        ! Local variables
        integer :: i, j, k
        real(dp) :: dif2, dself, dx, dy, dz, dk, aux, aux2

        !Mean-square displacement and intermediate scattering function
        dk = 6.6_dp
        do i = 1, nprom-1
            dif2 = 0.0_dp
            dself = 0.0_dp
            do j = 1, nprom-i
                do k = 1, np
                    dx = cfx(j+i,k)-cfx(j,k)
                    dy = cfy(j+i,k)-cfy(j,k)
                    dz = cfz(j+i,k)-cfz(j,k)
                    aux = dx*dx + dy*dy + dz*dz
                    dif2 = dif2 + aux
                    aux2 = dk*dsqrt(aux)
                    aux = dsin(aux2)/aux2
                    dself = dself+aux
                end do
            end do
            aux2 = real(np*(nprom-i))
            dif2 = dif2/aux2
            dself = dself/aux2
            wt(i) = wt(i)+dif2
            ft(i) = ft(i)+dself
        end do
    end subroutine difusion
end module observables