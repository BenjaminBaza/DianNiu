function particle_mesh_generator(method)
	global p
	global sol

	if method==1
		if sol.particle_mesh_distribution_order==1
			dxsn        = p.R_s_n/(sol.part_nb_cell);
			dxsp        = p.R_s_p/(sol.part_nb_cell);
			sol.part_coord_n  = 0:dxsn:p.R_s_n;
			sol.part_coord_p  = 0:dxsp:p.R_s_p;
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
	else

		sol.part_coord_n  = zeros(1,sol.part_nb_cell+1);
		sol.part_coord_p  = zeros(1,sol.part_nb_cell+1);
		sol.adjusted_position_a = zeros(1,sol.part_nb_cell+1);

		for i=1:1:sol.part_nb_cell+1
			proportion=(i-1)/sol.part_nb_cell;
			if i==1
				sol.adjusted_position_a(i)=0;
			elseif i<=(sol.part_nb_cell+1)/2
				sol.adjusted_position_a(i)=power(proportion,sol.particle_mesh_distribution_order) *0.5/power(0.5,sol.particle_mesh_distribution_order);
			else
				sol.adjusted_position_a(i)=1-sol.adjusted_position_a(sol.part_nb_cell + 2 - i);
			end
			sol.part_coord_n(i) = sol.adjusted_position_a(i)*p.R_s_n;
			sol.part_coord_p(i) = sol.adjusted_position_a(i)*p.R_s_p;
		end

	end



end