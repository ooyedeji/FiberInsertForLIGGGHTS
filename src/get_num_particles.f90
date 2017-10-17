subroutine get_num_particles_bulk_density_target

	use global_props
	
	implicit none

	num_particles = ceiling(vt*bulk_density/(atom_density*vp))
	
end subroutine get_num_particles_bulk_density_target

subroutine get_num_particles_void_ratio_target

	use global_props
	
	implicit none
	
	num_particles = ceiling(vt/((void_ratio + 1.0)*vp))
	
end subroutine get_num_particles_void_ratio_target

subroutine get_num_particles_void_fraction_target

	use global_props
	
	implicit none

	num_particles = ceiling(vt*(1-void_fraction)/vp)
	
end subroutine get_num_particles_void_fraction_target
