subroutine get_cylinder_options

	use math
	use global_props
	use cylinder_insert_props
	
	implicit none

	character(100)::arg
	real::arg_val
	integer::i

	do i=1,7
		read(*,*)arg,arg_val
		select case (arg)
			case ('radius_in')
				radius_in = arg_val
			case ('radius_out')
				radius_out = arg_val
			case ('z_low')
				z_low = arg_val
			case ('z_high')
				z_high = arg_val
			case ('x_pos')
				x_pos = arg_val
			case ('y_pos')
				y_pos = arg_val
			case ('cylinder_direction')
				cylinder_direction = int(arg_val)
			case DEFAULT
				print*,'Bad geometry argument ',arg
				stop
		end select
	enddo
	
	vt = pi*(radius_out*radius_out - radius_in*radius_in)*(z_high-z_low)

end subroutine get_cylinder_options

subroutine get_box_options

	use global_props
	use box_insert_props
	
	implicit none
	
	character(100)::arg
	integer::i
	real::arg_val
	
	do i=1,7
		read(*,*)arg,arg_val
		select case (arg)
			case ('box_x_low')
				box_x_low = arg_val
			case ('box_x_high')
				box_x_high = arg_val
			case ('box_y_low')
				box_y_low = arg_val
			case ('box_y_high')
				box_y_high = arg_val
			case ('box_z_low')
				box_z_low = arg_val
			case ('box_z_high')
				box_z_high = arg_val
		end select
	enddo
	
	vt = (box_x_high-box_x_low)*(box_y_high-box_y_low)*(box_z_high-box_z_low)

end subroutine get_box_options
