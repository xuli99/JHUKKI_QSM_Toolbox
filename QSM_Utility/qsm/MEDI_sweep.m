function chi_res = MEDI_sweep(deltaB, Params, DPWeight, maskErode, lambda, merit, edgePer, outputFile)
% Do lambda sweep search for MEDI
% Calculate curverture and save intemedia results
% Using maximum curverture to determine the best lambda
    
    lambda = sort(lambda);  % change to assending order

    Nlambda = length(lambda);

    regv_f = zeros(Nlambda, 1);
    datav_f = zeros(Nlambda, 1);
    chi_res_sweep = repmat(zeros(size(deltaB)), [1,1,1,Nlambda,1]);
    
    for lambda_ii = 1:Nlambda
        [chi_res_sweep(:,:,:,lambda_ii,1), regv, datav] = ...
            delta2chi_MEDI(deltaB, Params, DPWeight, maskErode, lambda(lambda_ii), merit, edgePer);  
        index_f = find(regv > 0, 1, 'last');
        regv_f(lambda_ii) = regv(index_f);
        datav_f(lambda_ii) = datav(index_f);
        saveNII(squeeze(chi_res_sweep(:,:,:,lambda_ii,1)), [outputFile, '_lambda_', num2str(lambda(lambda_ii))], Params, 1);
    end

    % select best lambda
    [ Kappa ] = calc_curv_spline( lambda, regv_f, datav_f, false, 1);
    [~] = draw_lcurve_median( lambda, regv_f, datav_f, 100);
    
    lambda_index_opt = (Kappa == max(Kappa));
    chi_res = chi_res_sweep(:,:,:,lambda_index_opt, 1);
