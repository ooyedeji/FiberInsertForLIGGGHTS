subroutine write_liggghts(num_inserted)

	use math
	use global_props
	
	implicit none
	
	integer::num_inserted,i,j
	integer::nbond_type,ii,atom_id
	real::ave_mass
	
	open(unit=2,file="liggghts_input_deck.txt",status="unknown",action="write")
	
	write(2,*)"Generated File for use with LIGGGHTS gran/bond hybrid"
	write(2,*)num_inserted*num_atoms_per_particle,"atoms"
	write(2,*)num_atom_styles,"atom types"
	write(2,*)(num_atoms_per_particle-1)*num_inserted,"bonds"
	if(num_atoms_per_particle>1)then
		nbond_type = 1
	else
		nbond_type = 0
	endif
	write(2,*)nbond_type,"bond types"
	
	write(2,*)' '
	write(2,*)sim_box_x_low,sim_box_x_high,'xlo xhi'
	write(2,*)sim_box_y_low,sim_box_y_high,'ylo yhi'
	write(2,*)sim_box_z_low,sim_box_z_high,'zlo zhi'
	
	ave_mass = atom_density*pi*4./3.*particle_radius*particle_radius*particle_radius
	write(2,*)' '
	write(2,*)'Masses'
	write(2,*)' '
	do i = 1,num_atom_styles
		write(2,*)i,ave_mass
	enddo
	
	write(2,*)' '
	write(2,*)'Atoms'
	write(2,*)' '
	ii = 0
	do i = 1,num_inserted
		do j = 1,num_atoms_per_particle
			ii = ii + 1
			write(2,"(I9,I2,5ES28.20,I2)")ii,1,x(i,j),y(i,j),z(i,j),2.0*particle_radius,atom_density,1
		enddo
	enddo
	
	if(num_atoms_per_particle>1)then
		write(2,*)' '
		write(2,*)'Bonds'
		write(2,*)' '
		ii = 0
		atom_id = 0
		do i = 1,num_inserted
			do j = 1,num_atoms_per_particle
				atom_id = atom_id + 1
				if(j<num_atoms_per_particle)then
					ii = ii + 1
					write(2,"(I9,I2,2I9)")ii,1,atom_id,atom_id+1
				endif
			enddo
		enddo
	endif
	
	close(2)

end subroutine write_liggghts