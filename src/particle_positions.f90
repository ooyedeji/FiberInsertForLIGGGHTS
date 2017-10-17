subroutine particle_positions(particles_inserted)

	use global_props
	
	implicit none
	
	integer::particles_inserted
	real::vs
	
	! Allocate Arrays
	allocate( 	x(num_particles,num_atoms_per_particle),&
				y(num_particles,num_atoms_per_particle),&
				z(num_particles,num_atoms_per_particle))
	
	print*,'Trying to insert',num_particles,'particles'
	
	select case (insert_type)
		case (1)
			call particle_positions_cylinder(particles_inserted)
	endselect
	
	! Gather statistics about the insertion
	vs = particles_inserted*vp
	void_ratio = vt/vs - 1.0
	bulk_density = atom_density*vs/vt
	void_fraction = 1.0-vs/vt
	
	print*,' Bulk Density = ',bulk_density
	print*,'   Void Ratio = ',void_ratio
	print*,'Void Fraction = ',void_fraction
	


end subroutine particle_positions

subroutine particle_positions_cylinder(num_inserted)

	use math
	use global_props
	use cylinder_insert_props
	
	implicit none
	
	integer,parameter::max_failed = 25
	integer::cur_failed,particles_inserted,tries
	integer::i,ii,part,num_inserted

	real::max_particle_length,tol,z_loc,r_loc,ang,xy,xy_mid,xy2
	real::sep_vec(num_atoms_per_particle),x1(3),x2(3),x3(3),x4(3),rand_num(3)
	real::temp_1(2),temp_2(2),a1(3),a2(3),temp_0(2),temp_vec
!	real::void_ratio,bulk_density,vs

	logical::in_box,is_touching,skipped,test1,test2,test3,test4,test5,test6,test7
	
	cur_failed = 0
	particles_inserted = 0
	! max_particle_length = 2.0*particle_radius+real(num_atoms_per_particle-1)*bond_length
	max_particle_length = num_atoms_per_particle*2.0*particle_radius
	tol = 1.00000*(2.0*particle_radius)
	
	sep_vec(1) = particle_radius
	do i = 2,num_atoms_per_particle
		sep_vec(i) = sep_vec(i-1) + 2.0*bond_radius
	enddo
	
	do part = 1,num_particles
		do ii = 1,num_atoms_per_particle
			x(part,ii) = 0.0
			y(part,ii) = 0.0
			z(part,ii) = 0.0
		enddo
	enddo

	ii = 0
	part = 0
	in_box = .false.
	is_touching = .false.
	skipped = .false.
	
	do while (ii < num_particles)
		part = part + 1
		do tries = 1,max_insert_tries
	
			in_box = .false.
			is_touching = .false.
			skipped = .false.
			
			call random_number(rand_num)

			z_loc = (z_high-z_low-tol)*rand_num(1) + z_low + tol/2.0
			if(radius_in.gt.0.0)then
				r_loc = (radius_out-radius_in-tol)*rand_num(2) + radius_in + tol/2.0
			else
				r_loc = (radius_out-0.5*tol)*rand_num(2)
			endif
			ang = 2.*pi*rand_num(3)
			
			x1(1) = r_loc*cos(ang) + x_pos
			x1(2) = r_loc*sin(ang) + y_pos
			x1(3) = z_loc
			
			call random_number(rand_num)
			rand_num = 2.0*rand_num - 1.0
			rand_num = norm_vec(rand_num)
			x2(1) = x1(1) + max_particle_length*rand_num(1)
			x2(2) = x1(2) + max_particle_length*rand_num(2)
			x2(3) = x1(3) + max_particle_length*rand_num(3)
			
			temp_0(1) = x1(1)
			temp_0(2) = x1(2)
			temp_1(1) = x2(1)
			temp_1(2) = x2(2)
			temp_2(1) = x_pos
			temp_2(2) = y_pos
			
			xy = get_norm(temp_1-temp_2)
			xy2= get_norm(temp_0-temp_2)

			test1 = xy    < (radius_out-tol)
			test2 = xy2   < (radius_out-tol)
			test3 = x2(3) < (z_high-tol)
			test4 = x2(3) > (z_low+tol)
			
			if(radius_in>0.0)then
				call min_dist_to_center(x1,x2,x_pos,y_pos,xy_mid)
				test5 = xy     > (radius_in+tol)
				test6 = xy_mid > (radius_in+tol)
				test7 = xy2    > (radius_in+tol)
			else
				test5 = .true.
				test6 = .true.
				test7 = .true.
			endif
			
			if(test1) then
				if(test2) then
					if(test3) then	
						if(test4) then
							if(test5) then
								if(test6) then
									if(test7) then
										in_box = .true.
									endif
								endif
							endif
						endif
					endif
				endif
			endif
			
			if(in_box)then
				if(num_atoms_per_particle.gt.1)then
					do i=1,(ii-1)
						x3(1) = x(i,1)
						x3(2) = y(i,1)
						x3(3) = z(i,1)
						x4(1) = x(i,num_atoms_per_particle)
						x4(2) = y(i,num_atoms_per_particle)
						x4(3) = z(i,num_atoms_per_particle)
						call find_intersect(x1,x2,x3,x4,a1,a2)
						if(get_norm(a1-a2)<tol)then
							is_touching = .true.
							exit
						endif
					enddo
				else
					do i=1,(ii-1)
						x3(1) = x(i,1)
						x3(2) = y(i,1)
						x3(3) = z(i,1)
						if(get_norm(x1-x3)<tol)then
							is_touching = .true.
							exit
						endif
					enddo
				endif
				
				if(.not.is_touching)then
					exit
				endif
				
				if(tries >= max_insert_tries) then
					skipped = .true.
					cur_failed = cur_failed + 1
					print*,'Current Failed ==',cur_failed
					exit
				endif
				
			endif
			
		enddo

		if(.not.skipped.and.in_box) then
			ii = ii + 1
			if(modulo(ii,25).eq.0)then
				print*,ii,' out of ',num_particles ,' particles inserted' 
			endif
			do i=1,num_atoms_per_particle
				x(ii,i) = x1(1) + sep_vec(i)*rand_num(1)
				y(ii,i) = x1(2) + sep_vec(i)*rand_num(2)
				z(ii,i) = x1(3) + sep_vec(i)*rand_num(3)
				cur_failed = 0
			enddo
		endif
		
		if(cur_failed>max_failed)then
			exit
		endif
	
	enddo
	
	num_inserted = ii
	
	! Correct the orientation
	select case (cylinder_direction)
		case (1)
			! Rotate into x
			do part=1,num_inserted
				do ii=1,num_atoms_per_particle
					temp_vec = x(part,ii)
					x(part,ii) = z(part,ii)
					z(part,ii) = temp_vec
				enddo
			enddo
			
		case (2)
			! Rotate into y
			do part=1,num_inserted
				do ii=1,num_atoms_per_particle
					temp_vec = y(part,ii)
					y(part,ii) = z(part,ii)
					z(part,ii) = temp_vec
				enddo
			enddo
		case (3)
			! Were Done
	endselect
		
end subroutine particle_positions_cylinder

subroutine min_dist_to_center(x1,x2,x_pos,y_pos,min_dist)
	use math
	implicit none
	real::x1(3),x2(3),x_pos,y_pos
	real::min_dist
	real::px,py,n2,u,x,y,new_point(2)
	px = x2(1)-x1(1)
	py = x2(2)-x1(2)
	n2 = px*px + py*py
	u = ((x_pos-x1(1))*px + (y_pos-x1(2))*py) / n2
	u = max(min(u,1.0),0.0)
	x = x1(1) + u*px
	y = x1(2) + u*py
	new_point(1) = x-x_pos
	new_point(2) = y-y_pos
	min_dist = get_norm(new_point)
end subroutine min_dist_to_center

subroutine find_intersect(x1,x2,x3,x4,a1,a2)
	implicit none
	real,dimension(3)::x1,x2,x3,x4,a1,a2
	real,dimension(3)::x13,x43,x21
	real::d1343,d4321,d1321,d4343,d2121,denom,numer,mua,mub
	x13 = x1-x3
	x43 = x4-x3
	x21 = x2-x1
	d1343 = x13(1) * x43(1) + x13(2) * x43(2) + x13(3) * x43(3)
    d4321 = x43(1) * x21(1) + x43(2) * x21(2) + x43(3) * x21(3)
    d1321 = x13(1) * x21(1) + x13(2) * x21(2) + x13(3) * x21(3)
    d4343 = x43(1) * x43(1) + x43(2) * x43(2) + x43(3) * x43(3)
    d2121 = x21(1) * x21(1) + x21(2) * x21(2) + x21(3) * x21(3)
    denom = d2121 * d4343 - d4321 * d4321
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343
    mua = max(min(mua,1.0),0.0)
    mub = max(min(mub,1.0),0.0)
    a1 = x1 + mua*x21
    a2 = x3 + mub*x43
end subroutine find_intersect
