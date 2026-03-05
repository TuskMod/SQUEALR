#Parameters is list of parameters

SetVarParms <- function(parameters){

    # get the parameters with more than one value
    variable_messy <- parameters[which(lapply(parameters, length)>1)]
    # filter out parameters that are supposed to have more than one value (or that have values connected to values of other parameters)
    variable_messy <- variable_messy[names(variable_messy) %in% c('out.opts', 'input', names(variable_messy)[grep('^B1', names(variable_messy))], 'ss', 'mort_val_test') == FALSE]
    # get all combinations
    temptab <- expand.grid(variable_messy)
    # build out table of parameters that are defined in sync, e.g. B1, density, ss with specific contact
    canonical.params <- data.frame(contact = rep(parameters$contact, each = length(parameters$density)),
                                    density = rep(parameters$density, length(parameters$contact)),
                                    ss = rep(parameters$ss, length(parameters$contact)))

    B1 <- unlist(lapply(parameters$contact, function(x){parameters[paste0('B1__',x)]}))
    rename.to.contact <- function(vars, in.name='contact'){
        # pulls out contact names from parameter values separated by __ and sticks them in a table
        varname <- deparse(substitute(vars))
        names(vars) <- unlist(lapply(strsplit(names(vars), '__'), function(x) x[[2]]))
        names(vars) <- unlist(lapply(strsplit(names(vars), '\\d'), function(x) x[[1]]))
        tab <- data.frame(names(vars), vars)
        names(tab) <- c(in.name, varname)
        return(tab)
    }

    B1_tab <- as.data.table(rename.to.contact(B1))
    B1_tab[, density := rep(parameters$density, length.out = nrow(B1_tab))]
    canonical.params <- merge(canonical.params, B1_tab, on=contact)
    # join user defined variables with canonical parameters
    common.columns <- names(canonical.params)[names(canonical.params) %in% names(temptab)]
    result <- inner_join(temptab, canonical.params, by=common.columns)
    # calculate B2 (relative to B1)
    result$B2 <- result$B1*parameters$B2_B1_factor

    F1 <- unlist(parameters[grep('^F1__', names(parameters))])
    F2_int <- unlist(parameters[grep('^F2_int', names(parameters))])
    F2_B <- unlist(parameters[grep('^F2_B', names(parameters))])
    F2i_int <- unlist(parameters[grep('^F2i_int', names(parameters))])
    F2i_B <- unlist(parameters[grep('^F2i_B', names(parameters))])
    infpd <- unlist(parameters[grep('^infpd', names(parameters))])
    incub <- unlist(parameters[grep('^incub', names(parameters))])

    F1s <- rename.to.contact(F1)
    F2_ints <- rename.to.contact(F2_int)
    F2_Bs <- rename.to.contact(F2_B)
    F2i_ints <- rename.to.contact(F2i_int)
    F2i_Bs <- rename.to.contact(F2i_B)
    infpds <- rename.to.contact(infpd, 'variant')
    incubs <- rename.to.contact(incub, 'variant')

    # connect all the tables together on their "contact" column -- https://stackoverflow.com/a/34393416
    shaperate_table <- Reduce(function(x1, x2) merge(x1, x2, by='contact'), list(F1s, F2_ints, F2_Bs, F2i_ints, F2i_Bs))
    variant_table <- Reduce(function(x1, x2) merge(x1, x2, by='variant'), list(infpds, incubs))

    result <- inner_join(result, shaperate_table, by='contact')
    result <- inner_join(result, variant_table, by='variant')

    return(result)
}
