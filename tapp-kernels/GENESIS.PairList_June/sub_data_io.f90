
subroutine sub_i1_data_write(iunit, iarray, n)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_i4_data_write(iunit, iarray, n)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_r4_data_write(iunit, array, n)
        integer :: iunit, n
        real(4) :: array(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(1pe15.6)') (array(i),i=1,n)
        return
end subroutine


subroutine sub_i1_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i, nr
        read(iunit,'(i10)') nr
        read(iunit,'(i10)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_i4_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i, nr
        read(iunit,'(i10)') nr
        read(iunit,'(i10)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_r4_data_read(iunit, array, n)
        integer :: iunit, n
        real(4) :: array(n)
        integer :: i, nr
        read(iunit,'(i10)') nr
        read(iunit,'(1pe15.6)') (array(i),i=1,nr)
        return
end subroutine


