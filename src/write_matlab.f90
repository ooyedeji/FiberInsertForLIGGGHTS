subroutine write_matlab(num_inserted)
	
	use global_props
	
	implicit none
	
	integer::num_inserted,i,j
	
	open(unit=1,file="MATLAB_Output.txt",status="unknown",action="write")
	
	do i = 1,num_inserted
		do j = 1,num_atoms_per_particle
			write(1,*)x(i,j),y(i,j),z(i,j)
		enddo
	enddo
	
	close(1)

end subroutine write_matlab