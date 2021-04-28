        subroutine sub_idata_write(iunit, n, iarray)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i
        write(iunit,'(i8)') n
        write(iunit,'(i8)') (iarray(i),i=1,n)
        return
        end subroutine

        subroutine sub_idata_read(iunit, n, iarray)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i, nr
        read(iunit,'(i8)') nr
        read(iunit,'(i8)') (iarray(i),i=1,nr)
        return
        end subroutine

