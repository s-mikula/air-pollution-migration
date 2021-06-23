# Load some common libraries
library(tidyverse)
library(sf)
#library(rgdal)
library(broom)
library(lfe)
library(Formula)
#library(pscl)
library(sandwich)
library(lmtest)
#library(formula.tools)
library(ggspatial)

add_stars <- function(x, latex = TRUE, strict = FALSE){
  y <- rep("",length(x))
  
  if(strict){
    y[x < 0.05]  <- "*"
    y[x < 0.01] <- "**"
    y[x < 0.001] <- "***" 
  }else{
    y[x < 0.1]  <- "*"
    y[x < 0.05] <- "**"
    y[x < 0.01] <- "***"   
  }
  
  if(latex){
    y <- str_c("^{",y,"}")
  }
  
  return(y)
}

print_latex <- function(x, path, booktabs = TRUE, pos = TRUE){
  
  
  if("stat" %in% names(x)){
  
  if(booktabs){
  x <- x %>% 
    mutate(
      row_num = row_number(),
      term = ifelse(stat == "estimate" & row_num != 1, str_c("\\addlinespace ",term),term),
      term = ifelse(stat == "adj.r.squared", str_c("\\midrule ",term),term),
    ) %>% 
    select(-row_num,-stat)
  }else{
    x <- x %>% select(-stat)
  }
  
  }
  
  if(pos){print(x)}
    
  write.table(x, 
              file = path, 
              quote = FALSE, 
              sep = "\t&\t", 
              eol = "\\\\\n", 
              col.names = FALSE, 
              row.names = FALSE,
              na = ""
              )
}

get_felm_summary <- function(x){
  fs <- summary(x, robust = TRUE)
  
  cfs <- fs$coefficients %>% as_tibble()
  names(cfs) <- c("est","se","tstat","pval")
  cfs <- cfs %>% 
    rowwise() %>% 
    mutate(
      est = format(est, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE),
      se = format(se, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
        str_c("(",.,")",add_stars(pval, latex = FALSE, strict = TRUE))
    ) %>% 
    ungroup() %>% 
    mutate(
      term = row.names(fs$coefficients)
    ) %>% 
    select(term, est, se) %>% 
    gather(stat,value,-term) %>% 
    arrange(term,stat)
  
  fe <- getfe(x) %>% 
    select(term = fe) %>%
    distinct() %>% 
    mutate(
      stat = "fe",
      value = "Yes"
    ) %>% 
    as_tibble() %>% 
    mutate_all(as.character)
    
  
  ms <- tribble(
    ~term, ~stat, ~value,
    "Observations", "n", format(x$N,big.mark = ","),
    "Adjusted R2", "r2", format(fs$adj.r.squared, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE)
  )
  
  out <- bind_rows(cfs,fe,ms)
  colnames(out)[3] <- as.character(x$lhs)
  
  return(out)
}

get_felm_table <- function(x, labs = NULL, 
                           olabs = c("Settlement FE" = "KOD_OBEC", 
                                     "Year FE" = "year",
                                     "Observations" = "Observations",
                                     "Adjusted R2" = "Adjusted R2")){
  x <- x %>% 
    reduce(left_join, by = c("term","stat"))
  
  x1 <- x %>% 
    filter(stat %in% c("est","se","se_conley"))
  
  if(!is.null(labs)){
    if(length(labs)==1 & all(labs == "yby")){
      labs <- create_labs(x)
    }
  }
  
  if(!is.null(labs)){
  x1 <- x1 %>% 
    mutate(
      term = factor(term, levels = labs, labels = names(labs))
    )
  }
  
  x1 <- x1 %>% 
    arrange(term,stat) %>% 
    mutate_all(as.character) %>% 
    mutate(
      #term = ifelse(str_detect(stat,"se"),"",term)
      term = case_when(
        stat == "se" ~ "Clustered SE",
        stat == "se_conley" ~ "Conley SE",
        TRUE ~ term
      )
    )
  
  x2 <- x %>% 
    filter(!(stat %in% c("est","se","se_conley")))
  
  if(!is.null(labs)){
    x2 <- x2 %>% 
      mutate(
        term = factor(term, levels = olabs, labels = names(olabs))
      )
  }
  
  x2 <- x2 %>% 
    arrange(term,stat) %>% 
    mutate_all(as.character)
  
  bind_rows(x1,x2) %>% 
    select(-stat)
}

create_labs <- function(x){
  raw <- x$term %>% unique()
  
  year <- str_extract(raw,"\\d{4}")
  so2 <- str_extract(raw,"avg\\d{2}")
  so2[str_detect(so2,"94")] <- "1994"
  so2[str_detect(so2,"00")] <- "2000"
  
  names(raw) <- str_c("SO2 level in ",so2," * (year = ",year,")")
  
  return(raw)
}

felmsp <- function(msp, df, spatial = TRUE){
  bs <- msp %>% felm(data = df)
  
  if(spatial){
  
  bscol <- felm_column(bs, spatial = TRUE)
  
  spse <- conleySE(msp, df) %>% 
    mutate(
      p.valueSA = add_stars(p.valueSA)
    ) %>% 
    rowwise() %>% 
    mutate(
      across(
        where(is.double),
        format, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE
        )
    ) %>% 
    ungroup() %>% 
    mutate(
      estimate = str_c("[",std.errorSA,"]",p.valueSA),
      stat = "std.error.conley"
    ) %>% 
    select(term,stat,estimate)
  
  names(spse)[3] <- bs$lhs
  
  spse %>% 
    filter(term %in% bscol$term) %>% 
    bind_rows(bscol,.) %>% 
    mutate(
      rows = case_when(
        stat %in% c("estimate","std.error","std.error.conley") ~ "a_coefs",
        stat == "fe" ~ "b_fe",
        stat == "adj.r.squared" ~ "c_r2",
        stat == "obs" ~ "d_obs"
      ),
      stat = factor(stat, levels = c("estimate","std.error","std.error.conley","fe","adj.r.squared","obs"))
    ) %>% 
    arrange(rows,term,stat) %>% 
    mutate(
      across(.cols = everything(), as.character)
    ) %>% 
    select(-rows)
  }else{
    felm_column(bs, spatial = FALSE)
  }
}

felm_column <- function(x, spatial = FALSE){
  
  if(!spatial){
  COEFS <- x %>% 
    tidy() %>% 
    rowwise() %>% 
    mutate(
      estimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
        str_c(.,add_stars(p.value, strict = FALSE)),
      std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
        str_c("(",.,")")
    ) %>%
    ungroup() %>%
    select(term,estimate,std.error) %>% 
    gather(stat,value,-term) %>% 
    arrange(term,stat) 
  }else{
    COEFS <- x %>% 
      tidy() %>% 
      rowwise() %>% 
      mutate(
        estimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE), 
        std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
          str_c("(",.,")") %>% 
          str_c(.,add_stars(p.value, strict = FALSE))
      ) %>%
      ungroup() %>%
      select(term,estimate,std.error) %>% 
      gather(stat,value,-term) %>% 
      arrange(term,stat) 
  }
  
  FES <- x$fe %>% names()
  FES <- tibble(
    term = FES,
    stat = rep("fe",length(FES)),
    value = rep("Yes",length(FES))
  )
  
  STATS <- x %>% 
    glance() %>% 
    select(adj.r.squared) %>% 
    gather(stat,value) %>% 
    mutate(
      value = format(value,digits=1,nsmall=3,trim=TRUE),
      term = stat
    )
  
  STATS <- tibble(
    term = "Observations",
    stat = "obs",
    value = x$N %>% as.integer() %>% format(big.mark=",")
  ) %>% 
    bind_rows(
      STATS,.
    )
  
  out <- bind_rows(
    COEFS,FES,STATS
  ) %>% 
    select(term,stat,value)
  
  names(out)[3] <- x$lhs
  
  return(out)
}

delete_rows <- function(x, fe = FALSE){
  x <- x %>% 
    mutate(
      term = ifelse(str_detect(term,"prior1"),str_replace(term,"prior1","prior"),term)
    ) %>% 
    filter(!str_detect(term,"year:periodprior:rdd_seg")) %>% 
    filter(!str_detect(term,"year:periodafter:rdd_seg")) %>% 
    filter(!str_detect(term,"year:periodprior:area_seg")) %>% 
    filter(!str_detect(term,"year:periodafter:area_seg")) %>% 
    mutate(
      term = case_when(
        str_detect(term,"so2avg94\\)") & stat == "estimate" ~ "Period (2000-2015)*SO2 concentration (1994)",
        term == "KOD_OBEC" ~ "Municipality FE",
        term == "year" ~ "Year FE",
        term == "YOKRES" ~ "District * Year FE",
        term == "adj.r.squared" ~ "Adjusted R2",
        term == "Observations" ~ "Observations",
        term == "so2" ~ "\\ce{SO2} concentration (\\so)",
        term == "so2avg9440" & stat == "estimate" ~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{40}{\\so}",
        term == "so2avg9450" & stat == "estimate" ~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{50}{\\so}",
        term == "so2avg9460" & stat == "estimate" ~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{60}{\\so}",
        #str_detect(term,"after") & str_detect(term,"benefits") ~ "Period (2000-2015)*benefits",
        str_detect(term,"so2avg94 == 30") & !str_detect(term,"benefits") & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{30}{\\so}",
        str_detect(term,"so2avg94 == 40") & !str_detect(term,"benefits") & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{40}{\\so}",
        str_detect(term,"so2avg94 == 50") & !str_detect(term,"benefits") & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{50}{\\so}",
        str_detect(term,"so2avg94 == 60") & !str_detect(term,"benefits") & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{60}{\\so}",
        str_detect(term,"so2avg94 == 30") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg94 == 40") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg94 == 50") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg94 == 60") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg94 >= 40") & !str_detect(term,"benefits") & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $\\geq$ \\SI{40}{\\so}",
        str_detect(term,"so2avg94 >= 40") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg94 == 30") & str_detect(term,"benefits")  & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{30}{\\so}",
        str_detect(term,"so2avg94 == 40") & str_detect(term,"benefits")  & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{40}{\\so}",
        str_detect(term,"so2avg94 >= 40") & str_detect(term,"benefits")  & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $\\geq$ \\SI{40}{\\so}",
        str_detect(term,"so2avg94 == 50") & str_detect(term,"benefits")  & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{50}{\\so}",
        str_detect(term,"so2avg94 == 60") & str_detect(term,"benefits")  & stat == "estimate"~ "Pre-desulfurization \\ce{SO2} concentration $=$ \\SI{60}{\\so}",
        str_detect(term,"so2avg94 == 30") & str_detect(term,"benefits")  & stat == "std.error"~ "\\quad\\quad\\quad $\\times$ Eligibility for benefits $\\times$ Post-desulfurization period",
        str_detect(term,"so2avg94 == 40") & str_detect(term,"benefits")  & stat == "std.error"~ "\\quad\\quad\\quad $\\times$ Eligibility for benefits $\\times$ Post-desulfurization period",
        str_detect(term,"so2avg94 >= 40") & str_detect(term,"benefits")  & stat == "std.error"~ "\\quad\\quad\\quad $\\times$ Eligibility for benefits $\\times$ Post-desulfurization period",
        str_detect(term,"so2avg94 == 50") & str_detect(term,"benefits")  & stat == "std.error"~ "\\quad\\quad\\quad $\\times$ Eligibility for benefits $\\times$ Post-desulfurization period",
        str_detect(term,"so2avg94 == 60") & str_detect(term,"benefits")  & stat == "std.error"~ "\\quad\\quad\\quad $\\times$ Eligibility for benefits $\\times$ Post-desulfurization period",
        term == "mean" ~ "Mean (migration rate)",
        term == "sd" ~ "Standard deviation (migration rate)",
        str_detect(term,"so2avg00 == 5") & !str_detect(term,"benefits") & stat == "estimate"~ "Post-desulfurization \\ce{SO2} concentration $=$ \\SI{5}{\\so}",
        str_detect(term,"so2avg00 == 10") & !str_detect(term,"benefits") & stat == "estimate"~ "Post-desulfurization \\ce{SO2} concentration $=$ \\SI{10}{\\so}",
        str_detect(term,"so2avg00 == 15") & !str_detect(term,"benefits")& stat == "estimate"~ "Post-desulfurization \\ce{SO2} concentration $=$ \\SI{15}{\\so}",
        str_detect(term,"so2avg00 == 20") & !str_detect(term,"benefits") & stat == "estimate"~ "Post-desulfurization \\ce{SO2} concentration $=$ \\SI{20}{\\so}",
        str_detect(term,"so2avg00 == 5") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg00 == 10") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg00 == 15") & !str_detect(term,"benefits")& stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"so2avg00 == 20") & !str_detect(term,"benefits") & stat == "std.error"~ "\\quad\\quad\\quad  $\\times$  Post-desulfurization period",
        str_detect(term,"edu_sec") & stat == "estimate"~ "Secondary educated (%)",
        str_detect(term, "edu_ter") & stat == "estimate"~ "Tertiary educated (%)",
        str_detect(term, "AGE") & stat == "estimate"~ "Age structure",
        str_detect(term, "Udos") & stat == "estimate"~ "Unemployment rate",
        TRUE ~ ""
      )
    ) #%>% 
    # mutate(
    #   term = ifelse(stat == "std.error","",term)
    # )
  
  if(!fe) x <- x %>% filter(stat != "fe")
  
  return(x)
}

tidyglance <- function(x, dig = 3){
  
  if(class(x) != "felm") stop("tidyglance() works only with objects of class felm.")
  
  require(tidyr)
  require(broom)
  require(dplyr)
  
  coefs <- tidy(x) %>% 
    mutate(
      p.value = add_stars(p.value)
    ) %>% 
    mutate_if(is.numeric,
              format,
              digits = 1,
              nsmall = dig,
              scientific = FALSE,
              trim = TRUE
    ) %>% 
    mutate(
      estimate = str_c(estimate,p.value),
      std.error = str_c("(",std.error,")")
    ) %>% 
    select(term,estimate,std.error) %>% 
    pivot_longer(-term, names_to = "stat", values_to = x$lhs)
  
  fe <- tibble(
    term = x$fe %>% names
  ) %>% 
    mutate(
      stat = "fe",
      cname = "Yes"
    ) %>% 
    rename_at(.vars = vars(cname), function(y) x$lhs)
  
  glance(x) %>% 
    mutate_if(is.numeric,
              format,
              digits = 1,
              nsmall = dig,
              scientific = FALSE,
              trim = TRUE
    ) %>% 
    mutate(
      Observations = format(x$N, big.mark = ",")
    ) %>% 
    select(adj.r.squared,Observations) %>% 
    mutate(
      stat = "rstats"
    ) %>% 
    pivot_longer(
      -stat,
      names_to = "term",
      values_to = x$lhs
    ) %>% 
    bind_rows(
      coefs,fe,.
    )
}

descstat <- function(vars, periods, xdata = rdata_mig, xformula = as.Formula(value ~ I(period == "after") | 0 | 0 | KOD_OBEC)){
  
  xformula_mean <- xformula %>% update(. ~ 1)
  
  bygroupdiff <- xdata %>% 
    nest_by({{ vars }}) %>%
    summarise(
      est = suppressWarnings(felm(xformula, data = data), classes = "error") %>% tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% 
    mutate(type = "diff")
  
  bygroup <- xdata %>% 
    nest_by({{ vars }},{{ periods }}) %>%
    summarise(
      est = suppressWarnings(felm(xformula_mean, data = data), classes = "error") %>% tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% 
    mutate(type = "mean")
  
  alldiff <- suppressWarnings(felm(xformula, data = xdata), classes = "error") %>% tidy() %>% slice_tail(n=1) %>% mutate(type = "diff")
  
  allgroup <- xdata %>% 
    nest_by({{ periods }}) %>%
    summarise(
      est = suppressWarnings(felm(xformula_mean, data = data), classes = "error") %>% tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% mutate(type = "mean")
  
  bind_rows(bygroup,bygroupdiff,alldiff,allgroup) %>% 
    mutate(
      p.value = add_stars(p.value)
    ) %>% 
    rowwise() %>% 
    mutate(
      across(where(is.double), ~format(.x, digits = 1, nsmall = 2, trim = TRUE, scientific = FALSE))
    ) %>% 
    ungroup() %>% 
    mutate(
      estimate = ifelse(type == "diff", str_c(estimate,p.value), estimate),
      std.error = str_c("(",std.error,")")
    ) %>% 
    select(-p.value,-statistic,-term)
}

descstat_lm <- function(vars, periods, xdata = rdata_mig, xformula = as.Formula(value ~ I(period == "after")), cls = TRUE){
  
  xformula_mean <- xformula %>% update(. ~ 1)
  
  bygroupdiff <- xdata %>% 
    nest_by({{ vars }}) %>%
    summarise(
      est = lm(xformula, data = data) %>% 
        coeftest(., vcov. = vcovCL(.,data$KOD_OBEC)) %>% 
        tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% 
    mutate(type = "diff")
  
  bygroup <- xdata %>% 
    nest_by({{ vars }},{{ periods }}) %>%
    summarise(
      est = lm(xformula_mean, data = data) %>% 
        coeftest(., vcov. = vcovCL(.,data$KOD_OBEC)) %>% 
        tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% 
    mutate(type = "mean")
  
  alldiff <- lm(xformula, data = xdata) %>% 
    coeftest(., vcov. = vcovCL(.,xdata$KOD_OBEC)) %>% 
    tidy() %>% slice_tail(n=1)  %>% mutate(type = "diff")
  
  allgroup <- xdata %>% 
    nest_by({{ periods }}) %>%
    summarise(
      est = lm(xformula_mean, data = data) %>% 
        coeftest(., vcov. = vcovCL(.,data$KOD_OBEC)) %>% 
        tidy() %>% slice_tail(n=1) %>% list(),
      .groups = "drop"
    ) %>% 
    unnest_wider(est) %>% mutate(type = "mean")
  
  bind_rows(bygroup,bygroupdiff,alldiff,allgroup) %>% 
    mutate(
      p.value = add_stars(p.value)
    ) %>% 
    rowwise() %>% 
    mutate(
      across(where(is.double), ~format(.x, digits = 1, nsmall = 2, trim = TRUE, scientific = FALSE))
    ) %>% 
    ungroup() %>% 
    mutate(
      estimate = ifelse(type == "diff", str_c(estimate,p.value), estimate),
      std.error = str_c("(",std.error,")")
    ) %>% 
    select(-p.value,-statistic,-term)
}

doublei <- function(x){
  
  y <- double(length(x))
  y[is.double(y)] <- NA 
  
  y[!str_detect(x,"i$")] <- as.double(x[!str_detect(x,"i$")])
  y[str_detect(x,"i$")] <- x[str_detect(x,"i$")] %>% 
    str_remove("[+-]?\\d+\\.\\d+i$") %>% 
    str_remove("[+-]?\\d+i$") %>% as.double()
  
  if(any(is.na(y))) warning("doublei() returns NAs.")
  
  return(y)
}

conleySE <- function(msp, dt){
  
  data_sp <- dt %>% 
    select(all_of(c(all.vars(msp),"refpoint_lon","refpoint_lat"))) %>% 
    drop_na()
  
  # Y
  Y <- data_sp %>% select(all_of(all.vars(msp)[1]))
  Y %>% write_csv(path = normalizePath("SpAdj/Ymat.csv"), col_names = FALSE)
  
  # Design matrix
  X <- model.matrix(msp, data = data_sp) %>% as_tibble()
  Xfe <- terms(msp, rhs = 2, lhs = 0) %>% as.character() %>% c(.,"-1") %>% str_c(collapse = "") %>% as.formula()
  Xfe <- model.matrix(Xfe, data = data_sp) %>% as_tibble()
  
  X <- bind_cols(X,Xfe)
  X %>% write_csv(path = normalizePath("SpAdj/Xmat.csv"), col_names = FALSE)
  
  data_sp %>% 
    select(refpoint_lat,refpoint_lon) %>% 
    write_csv(path = normalizePath("SpAdj/COORDmat.csv"), col_names = FALSE)
  
  data_sp %>% 
    select(KOD_OBEC) %>%
    mutate_all(as.integer) %>% 
    write_csv(path = normalizePath("SpAdj/IDmat.csv"), col_names = FALSE)
  
  data_sp %>% 
    select(year) %>%
    mutate_all(as.integer) %>% 
    write_csv(path = normalizePath("SpAdj/YEARmat.csv"), col_names = FALSE)
  
  routinepath <- normalizePath("SpAdj/routine.m")
  system(str_c("octave ",routinepath), ignore.stdout = TRUE)
  
  SA <- read_csv(normalizePath("SpAdj/OUTmat.csv"), col_names = FALSE)
  names(SA) <- c("estimateSA","std.errorSA","p.valueSA")
  
  tibble(term = names(X)) %>% 
    bind_cols(.,SA) %>% 
    mutate()
  
  SA %>% 
    mutate(across(
      where(is.character),
      doublei
    )) %>% 
    bind_cols(
      tibble(term = names(X)),
      .
    )
}