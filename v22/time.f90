program main
    implicit none
    integer :: dt(8)
    real    :: beats

    call date_and_time(values=dt)
    print '(i4, 5(a, i2.2))', dt(1), '/', dt(2), '/', dt(3), ' ', &
                              dt(5), ':', dt(6), ':', dt(7)

    beats = (dt(7) + ((dt(6) - dt(4) + 60) * 60) + (dt(5) * 3600)) / 86.4
    print '("Beats: @", f0.2)', beats
end program main
