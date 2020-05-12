function [optimized_fidelity, C_optimized, fidelity_per_sample, eta] = optimization(d, epsilon, eta, C1, N, p, batch_size, sample_size, state)

delta_f = 1; %Set initial value for fidelity difference
fidelity_per_sample = zeros(sample_size,1);
counter_1 = 1;

initial_fidelity = average_fidelity(d, C1, N, p, batch_size, state);
df_dc = finite_difference(d, epsilon, C1, N, p, batch_size, state);
C1 = permute(C1, [3 1 2]);
C1_update = permute((C1 + eta * df_dc), [3 2 1]);


while delta_f > 1E-5

    updated_fidelity = average_fidelity(d, C1_update, N, p, batch_size, state);
    
    fidelity_per_sample(counter_1) = updated_fidelity;
    
    counter_1 = counter_1 + 1;
    
    delta_f = (updated_fidelity - initial_fidelity);

    if delta_f > 0
        
        counter_2 = 1;
        C1_update_counter = zeros(sample_size, 1, length(C1)); 
        
        while delta_f > 0
            
            eta = 1.05 * eta;
            
            C1_update = permute((C1 + eta * df_dc), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C1_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = average_fidelity(d, C1_update, N, p, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    elseif delta_f < 0
        
        counter_2 = 1;
        C1_update_counter = zeros(sample_size, 1, length(C1));
        
        while delta_f < 0
            
            eta = 0.5 * eta;
            
            C1_update = permute((C1 + eta * df_dc), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C1_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = average_fidelity(d, C1_update, N, p, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    end
    
    C1_update = C1_update_counter((counter_2 - 1), :, :);
    
end

optimized_fidelity = average_fidelity(d, C1_update, N, p, batch_size, state);

C_optimized = permute(C1_update, [3 2 1]);
