function particle_mesh_generator()
	global p
	global sol

	if sol.particle_mesh_distribution_order==1
		sol.dxsn        = p.R_s_n/(sol.part_nb_cell);
		sol.dxsp        = p.R_s_p/(sol.part_nb_cell);
		sol.part_coord_n  = 0:sol.dxsn:p.R_s_n;
		sol.part_coord_p  = 0:sol.dxsp:p.R_s_p;
	else
		sol.part_coord_n  = zeros(1,sol.part_nb_cell+1);
		sol.part_coord_p  = zeros(1,sol.part_nb_cell+1);

		for i=1:1:sol.part_nb_cell+1
			proportion=(i-1)/sol.part_nb_cell;
			if i==1
				adjusted_position=0;
			else 
				adjusted_position=power(proportion,1/sol.particle_mesh_distribution_order);
			end
			sol.part_coord_n(i) = adjusted_position*p.R_s_n;
			sol.part_coord_p(i) = adjusted_position*p.R_s_p;
		end
	end



end