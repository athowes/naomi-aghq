#' Inputs for rho_t1_out
# vector<Type> rho_t1_out(plhiv_t1_out / population_t1_out);
# vector<Type> plhiv_t1_out(A_out * plhiv_t1)
# vector<Type> plhiv_t1(population_t1 * rho_t1);

#' Inputs for alpha_t1_out
# vector<Type> alpha_t1_out(artnum_t1_out / plhiv_t1_out);
# vector<Type> artnum_t1_out(A_out * artnum_t1);
# vector<Type> artnum_t1(population_t1 * prop_art_t1);

#' Inputs for lambda_t1_out
# vector<Type> lambda_t1_out(infections_t1_out / (population_t1_out - plhiv_t1_out));
# vector<Type> infections_t1_out(A_out * infections_t1);
# vector<Type> infections_t1(lambda_t1 * (population_t1 - plhiv_t1));

# Add to Table S2?
