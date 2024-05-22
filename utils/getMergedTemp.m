function T_f_out = getMergedTemp(mass_flow_rates, temperatures)
% determine the temp of the fluid stream that comes out after the mixing of
% a number of inlet streams
    
    assert(length(mass_flow_rates) == length(temperatures), 'Input arguments M and T must be of same length.')
    
    for i = 1:length(mass_flow_rates)
        h(i) = getSpecificEnthalpy(temperatures(i));
    end
    
    h_f_out = sum( h.*mass_flow_rates ) / sum( mass_flow_rates );
    T_f_out = getTemperature(h_f_out);
    
    
    %%%% subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function my_h = getSpecificEnthalpy(my_temp)
        
        % for water as compressed liquid, from themo textbook table A-4
        Tsat = [0.01 5 10 15 20];
        hf = [0.001 21.020 42.022 62.982 83.915];
        
        assert(my_temp < max(Tsat) && my_temp > min(Tsat), 'Your temp value is outside the range of currently programmed values.')
        my_h = interp1(Tsat, hf, my_temp);
    end

    function my_temp = getTemperature(my_h)
        
        % for water as compressed liquid, from themo textbook table A-4
        Tsat = [0.01 5 10 15 20];
        hf = [0.001 21.020 42.022 62.982 83.915];
        
        assert(my_h < max(hf) && my_h > min(hf), 'Your specific enthalpy value is outside the range of currently programmed values.')
        my_temp = interp1(hf, Tsat, my_h);
    end
    
end
