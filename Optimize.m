function [optimized_fidelity, C_optimized, fidelity_per_sample, eta] = Optimize(d, epsilon, eta, C1, C2, Permute, N, px, py, pz, batch_size, state)

delta_f = 1;
fidelity_per_sample = zeros(100,1);
counter_1 = 1;

initial_fidelity = Fidelity(d, C1, C2, Permute, N, px, py, pz, batch_size, state);
df_dc_1 = FiniteDifference(d, epsilon, C1, C2, Permute, N, px, py, pz, batch_size, state);
C1 = permute(C1, [3 1 2]);
C1_update = permute((C1 + eta * df_dc_1), [3 2 1]);

while delta_f > 1E-5
    
    updated_fidelity = Fidelity(d, C1_update, C2, Permute, N, px, py, pz, batch_size, state);
    
    fidelity_per_sample(counter_1) = updated_fidelity;
    
    counter_1 = counter_1 + 1;
    
    delta_f = (updated_fidelity - initial_fidelity);
    
    if delta_f > 0
        
        counter_2 = 1;
        C1_update_counter = zeros(100,1,length(C1));
        
        while delta_f > 0
            
            eta = 1.5 * eta;
            
            C1_update = permute((C1 + eta * df_dc_1), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C1_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = Fidelity(d, C1_update, C2, Permute, N, px, py, pz, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    elseif delta_f < 0
        
        counter_2 = 1;
        C1_update_counter = zeros(100,1,length(C1));
        
        while delta_f < 0
            
            eta = 0.5 * eta;
            
            C1_update = permute((C1 + eta * df_dc_1), [3 2 1]);
            
            C1_update_counter(counter_2, :, :) = C1_update;
            
            counter_2 = counter_2 + 1;
            
            updated_fidelity = Fidelity(d, C1_update, C2, Permute, N, px, py, pz, batch_size, state);
            
            delta_f = (updated_fidelity - initial_fidelity);
            
        end
        
    end
    
    C1_update = C1_update_counter((counter_2 - 1), :, :);
    
end

optimized_fidelity = Fidelity(d, C1_update, C2, Permute, N, px, py, pz, batch_size, state);

C_optimized = permute(C1_update, [3 2 1]);