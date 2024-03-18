#### MAASLIN2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# phyloseq is the input phyloseq to be processed.
# formula is the right side of the formula (the variables correlated with phyloseq abundances)
# abun_cut and prv_cut set Maaslin2's native abundance and prevalence filters.
# normalization and transform set Maaslin2's native data normaliztion/transformation settings.
# analysis_method is a native setting of Maaslin2 that dictates the analysis method used.
# (TSS, AST, and LM are recommended for microbiome data according to the tool authors.)
# max_significance specifies the significance cutoff.
# correction denotes the method used for p value correction.
# tax_level is the taxonomy level being analyzed.
tidy_maaslin = function(phyloseq,formula,
                        abun_cut = 0.0001, prv_cut = 0.1,
                        normalization = "TSS", transform = "AST", analysis_method = "LM",
                        max_significance = 0.05,correction = "BH",
                        tax_level){
  require(Maaslin2)
  
  # Glom to the desired level
  phyloseq = phyloseq %>% tax_glom(tax_level)
  
  # Adapt formula for maas input
  fixed_effects = formula %>% str_split('[+]') %>% lapply(str_trim) %>% unlist()
  
  # Remove NA values from variables of interest
  for(i in fixed_effects){
    phyloseq = eval(parse(text=paste0('phyloseq %>% subset_samples(!is.na(',i,'))')))
  }
  
  # Prep the abundances
  o = phyloseq@otu_table %>% as.data.frame()
  # Add X_ to front, will remove afterward. Otherwise Maaslin corrupts the names.
  rownames(o) = paste0('X_',rownames(o))
  
  # Prep the metadata 
  s = phyloseq %>% psmelt() %>% 
    dplyr::select(-c(OTU,Abundance,any_of(c('Kingdom','Phylum','Class','Order','Family','Genus','Species',colnames(phyloseq@tax_table@.Data))))) %>% 
    unique %>% `rownames<-`(.$Sample) %>% dplyr::select(-Sample)
  
  # Run Maaslin2
  set.seed(421)
  maas.sp = Maaslin2(
    input_data = o, 
    input_metadata = s, 
    output = 'temp', 
    min_abundance = abun_cut,
    min_prevalence = prv_cut,
    max_significance = max_significance,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    fixed_effects = fixed_effects,
    correction = correction,
    standardize = F,
    plot_heatmap = F,
    plot_scatter = F)
  
  # Removes Maaslin results folder, not needed
  unlink('temp',recursive = T) 
  
  # Select results of interest and format table
  ma = maas.sp$results %>% rename(OTU = feature,Variable = name, stat = coef, se = stderr,p_val = pval,q_val = qval) %>% 
    mutate(OTU = str_sub(OTU,start=3)) %>% # Remove 'X_' prefix
    # Add taxonomy info
    left_join(phyloseq@tax_table@.Data %>% as.data.frame %>% rownames_to_column('OTU') %>% 
                mutate(OTU = str_replace_all(OTU,'[^[:alnum:]|^_]','.'))) %>% 
    mutate(DA = 'Maaslin2') %>% dplyr::select(all_of(tax_level),DA,Variable,se,stat,p_val,q_val) %>% 
    pivot_longer(cols = c(se,stat,p_val,q_val),names_to = 'test',values_to = 'value') %>% 
    pivot_wider(names_from = c(Variable,test), names_sep = '__', values_from = value)
  
  return(list(results=ma,input_ps = phyloseq))
}

#### ANCOM-BC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# phyloseq is the input phyloseq to be processed.
# formula is the right side of the formula (the variables correlated with phyloseq abundances)
# p_adj_method denotes the method used for p value correction (fdr is equivalent to BH across these tools).
# prv_cut set ANCOM-BC's native prevalence filter.
# abun_cut sets the prevalence filter that will be applied after ANCOM-BC has been run (same as ALDEX2).
# group and struc_zero are native ANCOM-BC settings that are useful when the output variable is categorical (see documentation)
# tax_level is the taxonomy level being analyzed.
# alpha specifies the significance cutoff.
# normalization and transform set Maaslin2's native data normaliztion/transformation settings.
# analysis_method is a native setting of Maaslin2 that dictates the analysis method used.
# (TSS, AST, and LM are recommended for microbiome data according to the tool authors.)
tidy_ancombc = function(phyloseq, formula, p_adj_method = "fdr", prv_cut = 0.1, abun_cut = 1e-4,
                        group = NULL, struc_zero = FALSE,tax_level = NULL,
                        alpha = 0.05){
  
  require(ANCOMBC)
  
  # phyloseq=ps$rxn;formula='Status';group='Status';struc_zero=T; p_adj_method = "fdr"
  # abun_cut=1e-4;prv_cut=0.1; alpha = alpha
  # tax_level='Species'
  
  # Glom the pyloseq object
  phyloseq = phyloseq %>% tax_glom(tax_level)
  
  # Extract variables
  f = formula %>% str_split('[+]') %>% lapply(str_trim) %>% unlist()
  
  # Remove NA values from variables of interest
  for(i in f){
    phyloseq = eval(parse(text=paste0('phyloseq %>% subset_samples(!is.na(',i,'))')))
  }
  
  # Set seed and run ancom-bc
  set.seed(421); print('Starting ANCOM...')
  ancom.sp = ancombc(data=phyloseq,tax_level = NULL, formula = formula, p_adj_method = p_adj_method,
                     prv_cut = prv_cut,group = group, struc_zero = struc_zero, 
                     alpha = alpha)
  print('ANCOM done')
  
  # Abundance filter
  temp = phyloseq %>% microbiome::transform('compositional')
  #  Identify which taxa should be kept
  o = prune_taxa((taxa_sums(temp)/nsamples(temp))>abun_cut,temp) %>% psmelt() %>% .[[tax_level]] %>% unique
  an = ancom.sp$res %>% lapply(function(x)x %>% rownames_to_column('Taxon') %>% filter(Taxon %in% o))
  # Readjust p values
  an$q_val = an$p_val %>% mutate_at(vars(-contains('Taxon')),function(x)x %>% p.adjust(method='BH'))
  
  # ANCOM-BC output is a list of statistical tests - must be combined, then wrangled to match format of other tests
  for(i in 1:length(an)) an[[i]]$stat = names(an)[i]
  an = an %>% bind_rows %>% dplyr::select(-contains('Intercept')) %>%
    pivot_longer(cols = -c(Taxon,stat), names_to = 'Variable', values_to = 'value') %>% 
    mutate(stat = recode(.$stat, p = 'p_val',q = 'q_val',diff = 'diff_abn'),
           DA = 'ANCOM-BC') %>% 
    pivot_wider(names_from = c(Variable,stat), names_sep = '__', values_from = value)
  names(an)[1] = tax_level
  
  return(list(results = an,bias = ancom.sp$samp_frac, input_ps = phyloseq))
}

#### ALDEx2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

# phyloseq is the input phyloseq to be processed.
# formula is the right side of the formula (the variables correlated with phyloseq abundances)
# p_adj_method denotes the method used for p value correction (fdr is equivalent to BH across these tools).
# prv_cut and abun_cut sets the abundance and prevalence filters that will be applied after ALDEx2 has been run (same as ANCOM-BC abundance filter).
# tax_level is the taxonomy level being analyzed.
# max_significance specifies the significance cutoff.
tidy_aldex2 =  function(phyloseq,formula,prv_cut = 0.10,abun_cut=0.0001,
                        tax_level,max_significance){
  require(ALDEx2)

  # phyloseq=p;formula='Status';
  # abun_cut=1e-4;prv_cut=0.33; max_significance = alpha;
  # tax_level='Phylum'
  
  # glom the phyloseq object
  phyloseq = phyloseq %>% tax_glom(tax_level)
  
  # Get rid of integers, if necessary. Integers are found in functional data.
  # Integers can be safely rounded because they contribute 2% variance at most per functional annotation.
  # Calculation: the smallest nonzero mean functional abundance is ~45 reads. 45 can therefore be anything from 44.5-45.49. 1/45 = 0.022, 2.2% variance. Furthermore, these low abundance functions are filtered from the results, so the effective influence on variance is much lower than 2.2%.
  phyloseq = transform_sample_counts(phyloseq, function(x) round(x,0))
  
  ## Taxonomic level
  t = which(!is.na(phyloseq@tax_table[1,])) %>% max # select lowest taxonomic level
  t = colnames(phyloseq@tax_table)[t] # Name of tax level - should be the same as tax_level
  t2 = str_sub(t,1,1) %>% str_to_lower() # ex. 's' for Species
  
  # Extract variables so that they can be identified
  f = formula %>% str_split('[+]') %>% lapply(str_trim) %>% unlist()
  
  # Remove NA values from variables of interest
  for(i in f){
    phyloseq = eval(parse(text=paste0('phyloseq %>% subset_samples(!is.na(',i,'))')))
  }
  
  # Sample data
  s = phyloseq %>% psmelt() %>% 
    dplyr::select(Sample,all_of(f)) %>% 
    unique
  if('depth' %in% names(s)){
    s = s %>% mutate(depth = scale(s$depth)[,1])
  }
  
  # reorder s to match phyloseq
  temp = tibble(Sample = sample_names(phyloseq),n = 1:nsamples(phyloseq))
  s = s %>% left_join(temp) %>% arrange(n) %>% 
    mutate(Sample = str_replace_all(Sample,'-','.')) %>% 
    column_to_rownames('Sample')
  
  # Create the model matrix
  if(formula=='Status'){ 
    m = s[[f]] 
  } else { 
    m=eval(parse(text=paste("model.matrix(~",formula,",data=s)",sep='')))
  }

  # Run ALDEX2
  if(formula=='Status') testtype = 't' else testtype = 'glm'
  # Perform ALDEx2 analysis
  set.seed(421)
  aldex.sp = aldex(data.frame(otu_table(phyloseq)), 
             conditions = m,
             test=testtype,
             effect=TRUE,
             include.sample.summary=F,
             denom="all") %>% 
    rownames_to_column('OTU') %>% 
    left_join(data.frame(tax_table(phyloseq)) %>% 
                dplyr::select(all_of(tax_level)) %>% rownames_to_column("OTU")) %>% 
    dplyr::select(OTU, all_of(tax_level),everything())
  
  if(testtype=='t'){
    aldex.sp = aldex.sp %>% dplyr::select(-we.ep,-we.eBH,-contains('rab.'),-contains('diff')) %>% 
      rename(p_val = wi.ep, q_val = wi.eBH,stat = effect,se=overlap) %>%
      pivot_longer(cols = -c(OTU,contains(tax_level)), values_to = 'value', names_to = 'stat') %>% 
      mutate(Variable = paste0(f,levels(m)[2]))
  } else {
    aldex.sp = aldex.sp %>%
      pivot_longer(cols = contains('model'), values_to = 'value', names_to = 'Variable') %>% 
      mutate(stat = ifelse(str_detect(Variable,'Estimate'),'stat',
                           ifelse(str_detect(Variable,'Std..Error'),'se',
                                  ifelse(str_detect(Variable,'t.value'),'t value',
                                         ifelse(str_detect(Variable,'Pr...t...BH'),'q_val',
                                                ifelse(str_detect(Variable,'Pr...t..'),'p_val',
                                                       NA)))))) %>% 
      mutate(Variable = str_remove_all(Variable,paste0(c('model.','.Estimate','.Std..Error','.t.value','.Pr...t...BH','.Pr...t..'),collapse = '|'))) %>% 
      dplyr::select(-contains('Intercept'))
  }

  # Arrange stats
  aldex.all = aldex.sp %>% 
    pivot_wider(names_from = stat, values_from = value) %>% 
    mutate(DA = 'ALDEx2', .after=tax_level)
  
  # Abundance filter - identify taxa to keep
  p2 = phyloseq %>% microbiome::transform('compositional')
  a = prune_taxa((taxa_sums(p2)/nsamples(p2))>abun_cut,p2) %>% taxa_names
  # Prevalence filter - identify taxa to keep
  p = apply(phyloseq@otu_table@.Data, 1, function(x) (sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))>prv_cut)
  p2 = rownames(phyloseq@tax_table@.Data)[p]
  # Filter taxa
  aldex.all.final = aldex.all %>% 
    filter(OTU %in% a & OTU %in% p2) %>% 
    group_by(Variable) %>% mutate(q_val = p.adjust(p_val, method='fdr')) %>% ungroup %>% 
    dplyr::select(-OTU) %>% 
    # Format to match other tools
    pivot_longer(cols = -c(contains(tax_level),DA,Variable),names_to = 'stat', values_to = 'value') %>% 
    pivot_wider(names_from = c(Variable,stat), values_from = value, names_sep = '__')
  
  return(list(results = aldex.all.final,input_ps = phyloseq))
}

#### dplyr::select SIG UNIVAR TAXA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# This function filters the results of each differential abundance tests to include only
# the taxa that were univariately significant.
# This is because our selection criteria for our taxa of interest require that the microbe is
# univariately and multivariately significant (unless otherwise noted).
# res = list of all the tool outputs.
# sig = vector of the significant taxa
select_univar_taxa = function(res,sig){
  n = names(res[[1]]$results[1]) # Name of first column (depends on tax level)
  L = res %>% 
    lapply(function(x){
      x$results = x$results %>% filter(!!sym(n) %in% sig) %>% 
        pivot_longer(cols = contains('__'),names_sep = '__',
                     names_to = c('Variable','stat'),values_to = 'value')
      stats = unique(x$results$stat)
      x$results = x$results %>% 
        pivot_wider(names_from = stat, values_from = value) %>% 
        group_by(DA, Variable) %>% mutate(q_val = p.adjust(p_val,method='BH')) %>% ungroup %>%
        pivot_longer(cols = stats, names_to = 'stats', values_to = 'value') %>% 
        pivot_wider(names_from = c(Variable,stats), names_sep = '__', values_from = value)
      return(x)})
  return(L)
}

#### COMBINE RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# res is a list of all the tool outputs.
# mintools is the minimum number of tools that must be significant for the microbe to be considered significant.
# alpha is the significance cutoff.
# stat denotes whether the significance cutoff applies to the p value or the q value.
combine_da_results = function(res,mintools = 2,alpha=0.05,stat='q_val'){
  
  # res=L;mintools=2;alpha=alpha;stat='q_val'
  
  # Use the first of the results tables in res to start building the results table
  all = res[[1]]$results
  for(i in 2:length(res)) all = full_join(all,res[[i]]$results) 
  
  # Filter to only include features that are present across all tool
  # (in case of any differences in application of the abundance and prevalence filters)
  n = all[,1] %>% table() %>% as.data.frame %>% 
    filter(Freq==3) # must be in all 3 tools
  all = all[sapply(as.vector(all[,1]),function(x) x %in% as.vector(n[,1])),]
  
  # Identify the association direction of each result
  all2 = all %>% dplyr::select(-contains('diff_abn'),-contains('.high'),-contains('.low')) %>% 
    pivot_longer(cols = contains(c('p_val','q_val','lfc','beta','stat','effect')),
                 names_to = c('Variable','Statistic'),
                 values_to = 'Value',
                 names_sep = '__') %>% 
    mutate(Statistic = replace(Statistic,
                                Statistic %in% c('stat','beta','effect','lfc'),'association_dir')) %>%
    filter(!is.na(Value) & Statistic == 'association_dir') %>% 
    mutate(Value = replace(Value, Value>0,'Positive'),
           Value = replace(Value, Value<0,'Negative')) %>% 
    dplyr::select(-contains('__')) %>% 
    unique() %>% dplyr::select(-Statistic) %>% 
    rename(association_dir = Value)
  
  # Wrangle results, add association direction
  all = all %>% dplyr::select(-contains('__W'),-contains('diff_abn')) %>% 
    pivot_longer(cols = contains(c('p_val','q_val','lfc','se','beta','stat','stat.low','stat.high','effect')),
                 names_to = c('Variable','Statistic'),
                 values_to = 'Value',
                 names_sep = '__') %>% 
    mutate(Statistic = replace(Statistic,
                               Statistic %in% c('stat','beta','effect','lfc'),'estimate')) %>%
    # aldex.kw() doesn't give a regular standard deviation like the others. Instead, it gives confidence intervals via aldex.effect().
    # These are represented by estimate.high and low.
    mutate(Statistic = replace(Statistic, Statistic == 'stat.low','estimate.low')) %>%
    mutate(Statistic = replace(Statistic, Statistic == 'stat.high','estimate.high')) %>%
    filter(!is.na(Value)) %>% 
    pivot_wider(names_from = Statistic, values_from = Value) %>% 
    suppressWarnings() %>% suppressMessages() %>% 
    left_join(all2)
  
  # Determine which microbes are significant, and in how many datasets
  if(stat=='q_val'){
    sig = all %>% mutate(sig = as.numeric(q_val<alpha))
    all2 = all %>% filter(q_val<alpha)
  } else if(stat=='p_val'){
    sig = all  %>% mutate(sig = as.numeric(p_val<alpha))
    all2 = all %>% filter(p_val<alpha)
  } else {print('stat must be "q_val" or "p_val".')}
  
  sig = sig %>% 
    dplyr::select(all_of(names(all)[1]), DA,Variable,sig) %>%
    pivot_wider(names_from = DA, values_from = sig)
  # Total # of significant tools per microbe
  sig$total_sig = apply(sig[,3:ncol(sig)],1,sum)
  # sig_in_at_least_n tells us if the microbe is significant in at least n tools.
  sig = sig %>% mutate(sig_in_at_least_n = as.numeric(total_sig>=mintools))
  
  # List of all the sig microbes
  sig_microbes = sig %>% filter(!is.na(sig_in_at_least_n)) %>% 
    .[.$sig_in_at_least_n==1,] %>% 
    left_join(all2 %>% dplyr::select(names(all2)[1],Variable,association_dir) %>% unique)
  
  return(list(combined_res = all %>% arrange(p_val) %>% arrange(q_val), 
              significance_summary = sig %>% arrange(-sig_in_at_least_n), 
              sig_hits = sig_microbes))
}
