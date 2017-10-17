subroutine read_in

	use math
	use global_props
	use output_options
	
	implicit none
	
	logical::eof
	character(100)::arg
	real::arg_value
	
	write_to_matlab = .false.
	write_to_liggghts = .false.
	rand_num_seed = 1.0
	eof = .false.
	do while(.not.eof)
		read(*,*)arg,arg_value
	
		select case (arg)
			
			case ('eof')
				eof = .true.
				
			case ('rand_num_seed')
				rand_num_seed(1) = int(arg_value)
			
			case ('write_to_matlab')
				if(arg_value.gt.0.0)then
					write_to_matlab = .true. 
				else
					write_to_matlab = .false.
				endif
			
			case ('write_to_liggghts')
				if(arg_value.gt.0.0)then
					write_to_liggghts = .true. 
				else
					write_to_liggghts = .false.
				endif
			
			case ('num_atom_styles')
				num_atom_styles = int(arg_value)
			
			case ('num_atoms_per_particle')
				num_atoms_per_particle = int(arg_value)
			
			case ('num_particles')
				num_particles = int(arg_value)
			
			case ('max_insert_tries')
				max_insert_tries = int(arg_value)
			
			case ('insert_type')
				insert_type = int(arg_value)
				if(insert_type.eq.1) then
					! insert_type = 1
					call get_cylinder_options()
				elseif(insert_type.eq.2) then
					! insert_type = 2
					call get_box_options()
				endif
			
			case ('void_ratio')
				void_ratio = arg_value
			
			case ('bulk_density')
				bulk_density = arg_value
				
			case ('void_fraction')
				void_fraction = arg_value
			
			case ('particle_radius')
				particle_radius = arg_value
			
			case ('bond_radius')
				bond_radius = arg_value
				
			case ('bond_length')
				bond_length = arg_value
				
			case ('atom_density')
				atom_density = arg_value
				
			case ('sim_box_x_high')
				sim_box_x_high = arg_value
			
			case ('sim_box_x_low')
				sim_box_x_low = arg_value
				
			case ('sim_box_y_high')
				sim_box_y_high = arg_value
				
			case ('sim_box_y_low')
				sim_box_y_low = arg_value
				
			case ('sim_box_z_high')
				sim_box_z_high = arg_value
				
			case ('sim_box_z_low')
				sim_box_z_low = arg_value
				
			case DEFAULT
				print*,arg
				print*,'Is not a valid options'
				stop
			
		end select
	
	enddo
	
	vp = num_atoms_per_particle*4./3.*pi*particle_radius*particle_radius*particle_radius
	call check_input_data
	
	! Set random number generator
	call set_random_seed(rand_num_seed)

end subroutine read_in

subroutine set_random_seed(rand_num)

	integer:: i,j,rand_num
	integer, allocatable :: seeds(:)
	call random_seed(size = i)
	allocate(seeds(i))
	do j = 1,i
		seeds(j) = rand_num+j
	enddo
	call random_seed(put=seeds)
	deallocate(seeds)
	
end


