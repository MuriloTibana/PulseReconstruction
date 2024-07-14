% Return the basis vectors and corresponding energies for the initial and
% final states of Hydrogen
function [states,energies] = hydrogenBasis(r,N_states,l)

%============= Grid Parameters ==============
dr = r(2)-r(1);
temp = ones(length(r),1);

%=========== State Calculation ==============
% l = 1;
H = -0.5*(diag(temp(1:end-1),-1) + diag(-2*temp,0) + ...
    diag(temp(1:end-1),1))/dr.^2;
H = H + 0.5*diag(temp*l*(l+1)./(r'.^2));
H = H + diag(-1./r);
[states,energies] = eigs(H,N_states,'smallestreal', ...
    'Tolerance',1e-10,'MaxIterations',1000);

%================ Normalize =================
for state=1:N_states
    total = sum(abs(states(:,state)).^2)*dr;
    states(:,state) = states(:,state)./sqrt(total);
end % Loop over states

end % Function end