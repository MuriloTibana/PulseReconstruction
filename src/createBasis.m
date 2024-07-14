function createBasis(grid,directory)

%===== L=0 Block ============================
l = 0; N_states = 1;
[initial_state,initial_energy] = basisCalculation(grid,N_states,l);
save([directory '/l0_state.mat'],'initial_state','-mat');
save([directory '/l0_energy.mat'],'initial_energy','-mat');

%===== L=0 Block ============================
l = 1; N_states = grid(end);
[middle_states,middle_energies] = basisCalculation(grid,N_states,l);
save([directory '/l1_states.mat'],'middle_states','-mat');
save([directory '/l1_energies.mat'],'middle_energies','-mat');

%===== L=2 Block ============================
l = 2; N_states = grid(end);
[final_states,final_energies] = basisCalculation(grid,N_states,l);
save([directory '/l2_states.mat'],'final_states','-mat');
save([directory '/l2_energies.mat'],'final_energies','-mat');

%===== Save Grid Parameters =================
save([directory '/grid.mat'],'grid','-mat');

end % Function end