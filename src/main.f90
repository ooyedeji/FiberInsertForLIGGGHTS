Program main

	implicit none

	real::time_start,time_end
	integer::particles_inserted

	call cpu_time(time_start)
	
	! Read in options
	call read_in()
	
	! Find particle positions
	call particle_positions(particles_inserted)
	
	! Write data to file
	call write_output(particles_inserted)
	
	call cpu_time(time_end)
	print*,'Total packing time ==',(time_end-time_start)
	
end Program
