        subroutine sub_idata_write(iunit, n, iarray)
        integer :: iunit, n, iarray(n)
        integer :: i, nr
        write(iunit,'(i8)') nr
        !cx write(*,'(a,i8)') "reading data items:", nr
        write(iunit,'(i8)') (iarray(i),i=1,nr)
        return
        end subroutine

        subroutine sub_idata_read(iunit, n, iarray)
        integer :: iunit, n, iarray(n)
        integer :: i, nr
        read(iunit,'(i8)') nr
        !cx write(*,'(a,i8)') "reading data items:", nr
        read(iunit,'(i8)') (iarray(i),i=1,nr)
        return
        end subroutine

