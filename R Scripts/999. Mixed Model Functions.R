# 999. Mixed Model Functions.R

# This function runs the mixed model without controlling for any covariables.
# df is the input dataset.
# species is the species (baseline abundance) being correlated with disease severity/progression.
# variable is the disease severity metric of interest (ex. total UPDRS score)
run_lmer = function(df,species,variable){
 
  # df=pd;species=sig[1];variable='mds.total'
  
  # Make sure everyone has at least 2 visits
  temp = df %>% filter(Species==species) %>% group_by(Sample) %>% count() %>% filter(n<2)
  df = df %>% filter(!(Sample %in% temp$Sample))
  
  # Temp name for the variable so it can be referenced
  names(df)[which(names(df)==variable)] = 'test'
  # Set range to 0-1
  df$test = rescale(df$test)
  # Center around 0
  df$test = df$test-mean(df$test,na.rm=T)
  
  # Run mixed model with LMER package.
  # lmerTest must also be loaded so that stats can be calculated.
  # "months_last_visit|Sample" denotes that each sample will have its own slope and intercept (i.e. random effect).
  #     See https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf for more detail.
  a = lmer(formula = test ~ Abundance*months_last_visit + (months_last_visit|Sample), 
           data = df %>% filter(Species==species), # Test 1 species at a time
           REML=F, # See https://stats.stackexchange.com/questions/272633/how-to-decide-whether-to-set-reml-to-true-or-false/272654#272654
           # Extra necessary parameters, used recommended values
           control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) %>% 
    # Format results
    summary %>% .$coefficients %>% as.data.frame() %>% rownames_to_column('Variable') %>% 
    mutate(Taxon = species,Test=variable)
  
  # Change 'test' to something more descriptive
  a$Variable = sapply(a$Variable,function(x)str_replace_all(x,'test',variable))
  # allres has to be initialized before the function is called - results will be collected here
  allres <<- allres %>% rbind(a)
}

# This function is the same as run_lmer(), but includes covariables that are geared for differences in **disease progression**.
# (i.e. linear changes in disease severity scores over the course of the study)
# See run_lmer function for in-function commenting.
# df is the input dataset.
# species is the species (baseline abundance) being correlated with disease severity/progression.
# variable is the disease severity metric of interest (ex. total UPDRS score)
run_lmer_multi = function(df,species,variable){
  
  # Make sure everyone has at least 2 visits
  temp = df %>% filter(Species==species, !is.na(!!sym(variable))) %>% 
    group_by(Sample) %>% count() %>% filter(n<2) %>% ungroup
  df = df %>% filter(!(Sample %in% temp$Sample))
  
  baseline = df %>% arrange(redcap_event_name) %>% group_by(Sample) %>% 
    filter(redcap_event_name==redcap_event_name[1L]) %>% ungroup() %>% 
    select(Sample,all_of(variable)) %>% 
    `colnames<-`(c('Sample','baseline'))
  
  names(df)[which(names(df)==variable)] = 'test'
  
  df = df %>% left_join(baseline) #%>% 
    # mutate(test = test-baseline)
  
  df$test = rescale(df$test)
  df$baseline = rescale(df$baseline)
  df$test = df$test-mean(df$test,na.rm=T)
  df$baseline = df$baseline-mean(df$baseline,na.rm=T)
  
  # print(min(df$test,na.rm = T));print(max(df$test,na.rm = T))
  
  a = lmer(test ~ Sex + laxatives + depth + baseline + Abundance*months_last_visit + (months_last_visit|Sample), 
           data = df %>% filter(Species==species),
           REML=F, 
           control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) %>% 
    summary %>% .$coefficients %>% as.data.frame() %>% rownames_to_column('Variable') %>% mutate(Taxon = species,Test=variable)
  a$Variable = sapply(a$Variable,function(x)str_replace_all(x,'test',variable))
  allresmulti <<- allresmulti %>% rbind(a)
}


# This function changes the output statistical table to be more descriptive for plotting purposes.
# df is the stats table
# v is a list of the disease severity variables that you want to include in the final results.
edit_lmer = function(df,v){
  
  df = df %>% filter(Test %in% v)
  
  # Edit levels for plotting
  df$Test[df$Test=='mds.total'] = 'Total UPDRS'
  df$Test[df$Test=='mds1.total'] = 'UPDRS 1'
  df$Test[df$Test=='mds2.total'] = 'UPDRS 2'
  df$Test[df$Test=='mds3.total'] = 'UPDRS 3'
  df$Test[df$Test=='mds4.total'] = 'UPDRS 4'
  df$Test[df$Test=='levo.eq.dose'] = 'Levo. Eq. Dose'
  df$Test[df$Test=='bdi_total'] = 'Depression'
  df$Test[df$Test=='fss_total'] = 'Fatigue'
  df$Taxon = str_remove_all(df$Taxon,'s__') %>% str_replace_all('_',' ')
  
  df$Test = factor(df$Test,levels = c('Total UPDRS','UPDRS 1','UPDRS 2','UPDRS 3','UPDRS 4','Levo. Eq. Dose','Depression','Fatigue'))
  df = df %>% droplevels() %>% as.data.frame()
  
  return(df)
}

# This function plots the results of the mixed model analysis. Each species of interest is on the Y axis, and each disease severity metric is on the x axis. Each correlation is represented as one dot.
# The dot colour represents the effect size of the correlation.
# corrected P values are labelled on the dots using symbols.
# df is the stats table output from run_lmer().
# v is a vector of the variables to include in the figure.
# bsize is the base size of the plot.
plot_lmer = function(df,v,bsize=20){
  # df = allres.p; v=c('mds.total','mds1.total','mds2.total','mds3.total','mds4.total')
  
  # See function description above
  df = edit_lmer(df,v)
  
  # Arrange data
  temp = df %>% filter(Test==(levels(.$Test)[1])) %>% arrange(Estimate)
  df$Taxon = factor(df$Taxon, levels = temp$Taxon)
  
  # Plot
  p = df %>% 
    mutate(Qval = ifelse(qval>0.05,'',
                         ifelse(qval>0.01,'*',
                                ifelse(qval>0.001,'**',
                                       ifelse(qval<=0.001,'***',''))))) %>% 
    mutate(Qval2 = ifelse(qval>0.1,'',
                          ifelse(qval>0.05,'+',''))) %>% 
    ggplot(aes(Test,Taxon,col = Estimate)) +
    geom_point(size =14) +
    theme_classic(base_size = bsize) +
    scale_color_gradient2(low="darkblue", high="darkred", guide="colorbar") +
    geom_text(aes(label=Qval),size = 12, col = 'white',nudge_y = -0.25) +
    geom_text(aes(label=Qval2),size = 12, col = 'white',nudge_y = 0) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}
