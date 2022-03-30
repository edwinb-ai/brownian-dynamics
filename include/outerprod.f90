module outerprod
    use types

    implicit none
    public outerprod_d

contains

function outerprod_d(a, b)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: outerprod_d(size(a), size(b))

    outerprod_d = spread(a, dim=2, ncopies=size(b)) * &
                    spread(b, dim=1, ncopies=size(a))

end function outerprod_d

end module outerprod