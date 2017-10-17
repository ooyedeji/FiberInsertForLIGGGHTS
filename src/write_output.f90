subroutine write_output(num_inserted)

	use math
	use global_props
	use output_options
	
	implicit none
	
	integer::num_inserted
	
	if(write_to_matlab)then
		print*,' '
		print*,'Writting MATLAB file'
		call write_matlab(num_inserted)
	endif
	if(write_to_liggghts)then
		print*,' '
		print*,'Writting LIGGGHTS file'
		call write_liggghts(num_inserted)
	endif	
	
	print*,' '
	print*,num_inserted,' particles inserted'
	
	deallocate(x,y,z)
		
end subroutine write_output
