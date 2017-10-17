module global_props

	integer,dimension(1)::rand_num_seed

	integer::num_atom_styles
	integer::num_atoms_per_particle
	integer::num_particles
	integer::max_insert_tries
	integer::insert_type
	
	real::particle_radius
	real::bond_radius
	real::bond_length
	real::atom_density
	real::vp ! Volume of a single particle
	
	real::sim_box_x_low
	real::sim_box_y_low
	real::sim_box_z_low
	real::sim_box_x_high
	real::sim_box_y_high
	real::sim_box_z_high
	
	real:: vt ! Volume of the insertion geometry
	
	real::bulk_density,void_ratio,void_fraction
	
	real,dimension(:,:),allocatable::x,y,z
	
end module global_props

module cylinder_insert_props

	real::radius_in
	real::radius_out
	real::z_low
	real::z_high
	real::x_pos
	real::y_pos
	
	integer::cylinder_direction

end module cylinder_insert_props

module box_insert_props

	real::box_x_low
	real::box_y_low
	real::box_z_low
	
	real::box_x_high
	real::box_y_high
	real::box_z_high

end module box_insert_props

module output_options

	logical::write_to_matlab
	logical::write_to_liggghts
	logical::write_to_EDEM_h5

end module output_options

module math

	implicit none
	
	private
	
	public::get_norm,pi,norm_vec
	
	real,parameter::pi = 3.141592653589793
	
	contains 
	
		function get_norm(x)
		
			real::x(:),sum_x,get_norm
			integer::i,n
			
			n = size(x)
			sum_x = 0
			do i = 1,n
				sum_x = sum_x + x(i)*x(i)
			enddo
			get_norm = sqrt(sum_x)
			
		end function get_norm
		
		function norm_vec(x)
		
			real::x(:)

			real,dimension(1:size(x))::norm_vec
			
			norm_vec = x/get_norm(x)
			
		end function norm_vec
		
end module math
