function [optimized_fidelity, C_optimized, fidelity_per_sample, eta] = optimization(d, epsilon, eta, C_variable, C_constant, N, px, py, pz, batch_size, sample_size, state)

delta_f = 1;
fidelity_per_sample = zeros(sample_size,1);
counter_1 = 1;

initial_fidelity = average_fidelity(d, C_variable, C_constant, N, px, py, pz, batch_size, state);
df_dc_1 = finite_difference(d, epsilon, C_variable, C_constant, N, px, py, pz, batch_size, state);
C_variable = permute(C_variable, [3 1 2]);
C_variable_update = permute((C_variable + eta * df_dc_1), [3 2 1]);

while delta_f > 1E-5
    
    updated_fidelity = average_fidelity(d, C_variable_update, C_constant, N, px, py, pz, batch_size, state);
    
    fidelity_per_sample(counter_1) = updated_fidelity;
    
    counter_1 = counter_1 + 1;
    
    delta_f = (updated_fidelity - initial_fidelity);
    
    if delta_f > 0
        
        counter_2 = 1;
        C1_update_counter = zeros(sample_size,1,length(C_variable));
        
        while delta_f > 0
            
            eta = 1.05 * eta;
            
            C_variable_update = permute((C_variable + eta * df_dc_1), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C_variable_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = average_fidelity(d, C_variable_update, C_constant, N, px, py, pz, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    elseif delta_f < 0
        
        counter_2 = 1;
        C1_update_counter = zeros(sample_size,1,length(C_variable));
        
        while delta_f < 0
            
            eta = 0.5 * eta;
            
            C_variable_update = permute((C_variable + eta * df_dc_1), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C_variable_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = average_fidelity(d, C_variable_update, C_constant, N, px, py, pz, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    end
    
    C_variable_update = C1_update_counter((counter_2 - 1), :, :);
    
end

optimized_fidelity = average_fidelity(d, C_variable_update, C_constant, N, px, py, pz, batch_size, state);

C_optimized = permute(C_variable_update, [3 2 1]);