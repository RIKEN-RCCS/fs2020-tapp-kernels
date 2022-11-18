
subroutine sub_i1_data_write(iunit, iarray, n,a)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i
        character(20) :: a
        write(iunit,'(a20,i10)') a,n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_i4_data_write(iunit, iarray, n,a)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i
        character(20) :: a
        write(iunit,'(a20,i10)') a,n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_r4_data_write(iunit, array, n,a)
        integer :: iunit, n
        real(4) :: array(n)
        integer :: i
        character(20) :: a
        write(iunit,'(a20,i10)') a,n
        write(iunit,'(1pe25.15)') (array(i),i=1,n)
        return
end subroutine


subroutine sub_i1_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i, nr
        character(20) :: a
        read(iunit,'(a20,i10)') a,nr
        read(iunit,'(i10)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_i4_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i, nr
        character(20) :: a
        read(iunit,'(a20,i10)') a,nr
        read(iunit,'(i10)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_r4_data_read(iunit, array, n)
        integer :: iunit, n
        real(4) :: array(n)
        integer :: i, nr
        character(20) :: a
        read(iunit,'(a20,i10)') a,nr
        read(iunit,'(1pe25.15)') (array(i),i=1,nr)
        return
end subroutine


