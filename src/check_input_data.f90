subroutine check_input_data

	use global_props
	use output_options
	
	implicit none
	
	logical::end_program,has_output
	
	end_program = .false.
	
	if(num_atom_styles<1)then
		print*,'Number of atom styles must be at least 1'
		end_program = .true.
	endif
	
	if(num_atoms_per_particle<1)then
		print*,'Number of atoms per particle must be at least 1'
		end_program = .true.
	endif
	
	if(max_insert_tries < 10000)then
		if(max_insert_tries < 100)then
			print*,'Max insert tries must be greater than 10'
			end_program = .true.
		else
			print*,'WARNING: It is suggested to use more than 10000...'
		endif
	endif
	
	if(particle_radius<=0.0)then
		print*,'Your atom must have a radius...'
		end_program = .true.
	endif
	
	if((bond_radius<=0.0).and.(num_atoms_per_particle>1))then
		print*,'bond radius must be greater than 0'
		end_program = .true.
	endif
	
	if((bond_length<=0.0).and.(num_atoms_per_particle>1))then
		print*,'bond length must be greater than 0'
		end_program = .true.
	endif
	
	if(atom_density<=0.0)then
		print*,'Atoms must have a density greater than 0'
		end_program = .true.
	endif
	
	has_output = write_to_liggghts.or.write_to_matlab
	has_output = has_output.or.write_to_EDEM_h5
	if(.not.has_output)then
		print*,'Program must write output to file...'
		end_program = .true.
	endif
	
	if(num_particles.le.0)then
		if(void_ratio.gt.0.0)then
			call get_num_particles_void_ratio_target
		elseif(bulk_density.gt.0.0)then
			call get_num_particles_bulk_density_target
		elseif(void_fraction.gt.0.0)then
			call get_num_particles_void_fraction_target
		else
			print*,'You must set the number of particles or set a packing property'
			end_program = .true.
		endif
	endif
	
	if(end_program)then
		print*,'Stopping...'
		stop
	endif

end subroutine check_input_data
