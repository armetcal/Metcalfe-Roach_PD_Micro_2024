# Cross-sectional functions

# This function runs Spearman correlations between var and all significant taxa (one at a time),
# and combines the results into a table.
# df = the dataset to analyze
# var = the variable of interest
# grp = specify what subset of the total dataset is being tested (all, pd, ctrl, etc)
correlate_met = function(df,var,grp){
    # For every microbe of interest...
    test = apply(df[,sig],2,function(x){
      ct = cor.test(x,df[[var]],method='spearman') %>% broom::tidy() %>% as.data.frame() 
      # Extract CI from rank-transformed Pearson tests
      ci = cor.test(rank(x),rank(df[[var]]),method='pearson') %>% broom::tidy() %>% as.data.frame()
      ct = cbind(ct,ci[,colnames(ci) %in% c('conf.high','conf.low')]) %>% t() %>% as.data.frame()
      return(ct)
  }) %>% # Then format the results
      as.data.frame() %>% `colnames<-`(sig) %>% 
      t() %>% as.data.frame() %>% 
      `names<-`(c('Estimate','Statistic','Pval','Method','Alternative','Conf High','Conf Low')) %>% 
      rownames_to_column('taxon') %>% 
      # Add test info and correct p values
      mutate(variable = var,
             group = grp,
             Padj = p.adjust(.$Pval,method='fdr')) %>% 
      arrange(Pval) %>% 
      dplyr::select(group,variable,taxon,Estimate,`Conf Low`,`Conf High`,Pval,Padj,Statistic,everything())
  return(test)
}

# This function runs linear models between var and all significant taxa (one at a time),
# controlling for covariables, and combines the results into a table.
# It is very similar to the spearman function above.
# Difference: Rescaling species and metabolites to 0-1 before running so that all estimates are comparable!
#    This is not necessary for spearman, since rho values are not affected by the scale of the variable.
# df = the dataset to analyze
# var = the variable of interest
# grp = specify what subset of the total dataset is being tested (all, pd, ctrl, etc)
# Incl_lax is in reference to the fact that the metabolites do not correlate with laxative use independently of Bristol rating, but the metabolite~Bristol correlation is impacted by laxative use. If incl_lax is true, then laxative use will be included in the model as an interaction with bristol rating. If false, it will not be included as a variable. incl_lax=F should be used when laxative users are filtered out of the dataset.
glm_met = function(df,var,grp, incl_lax=T){

  # Create formula for glm
  if(grp=='All'){
    if(incl_lax==T){
      f = formula(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$Status+df$laxatives+df$bristol+df$depth)
    } else {
      f = formula(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$Status+df$bristol+df$depth)
    }
    
  } else {
    if(incl_lax==T){
      f = formula(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$laxatives+df$bristol+df$depth)
    } else {
      f = formula(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$bristol+df$depth)
    }
  }
  
  test = data.frame()
  for(i in 1:length(sig)){
    g = summary(glm(f))$coefficients
    g = g %>% as.data.frame %>% rownames_to_column('Explanatory')
    g$Explanatory[g$Explanatory=='rescale(df[[sig[i]]])'] = 'Microbe'
    g$Explanatory[g$Explanatory=='rescale(df[[sig[i]]]):df$StatusPD'] = 'Microbe*Status'
    test = test %>% rbind(g %>% mutate(Taxon = sig[i]))
  }
  test = test %>% 
    rename(pval = 'Pr(>|t|)') %>% 
    mutate(Explanatory = str_remove(.$Explanatory,'df$'),
           Response = var,group = grp) %>% 
    group_by(Explanatory) %>% 
    mutate(qval = p.adjust(pval,method = 'fdr')) %>% ungroup() %>% 
    arrange(pval) %>%
    select(group,Taxon,Response,Explanatory,everything())
  return(test)
}


glm_met_severity = function(df,var,grp){

  f = formula(rescale(df[[var]])~rescale(df[[sig[i]]]) + df$laxatives + df$depth + df$Sex + rescale(df$Disease.duration))

  test = data.frame()
  for(i in 1:length(sig)){
    g = summary(glm(f))$coefficients
    g = g %>% as.data.frame %>% rownames_to_column('Explanatory')
    g$Explanatory[g$Explanatory=='rescale(df[[sig[i]]])'] = 'Microbe'
    test = test %>% rbind(g %>% mutate(Taxon = sig[i]))
  }
  test = test %>%
    rename(pval = 'Pr(>|t|)') %>%
    mutate(Explanatory = str_remove(.$Explanatory,'df$'),
           Response = var,group = grp) %>%
    group_by(Explanatory) %>%
    mutate(qval = p.adjust(pval,method = 'fdr')) %>% ungroup() %>%
    arrange(pval) %>%
    select(group,Taxon,Response,Explanatory,everything())
  return(test)
}