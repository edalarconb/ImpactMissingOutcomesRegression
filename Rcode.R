
# ChJS ---------------------------------------------------------------
rm(list = ls())

# librerías ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,data.table,rio,openxlsx,kableExtra,gridExtra)


stoye_limits <- function(datos,Ymin,Ymax,cd,...){
  
  # * add column of 1 -------------------------------------------------------
  x = datos %>% 
    dplyr::mutate(const = 1) %>% 
    as_tibble()
  
  nt = dim(x)[1] #number of row data
  
  # * variable names ---------------------------------------------
  x.names.or <- names(x)
  names(x) =  letters[1:length(x.names.or)] #"a" will be the name of the response variable
  
  variables <- names(x)[-1]
  nx = length(variables) #number of X variables
  
  
  # * inputs for Ex under and upper -----------------------------------------
  y.min = Ymin
  y.max = Ymax
  
  
  # * unique x's ------------------------------------------------------------
  x.unique = x %>%
    dplyr::group_by_at(vars(variables)) %>% 
    dplyr::summarise(nd=n(),
                     ynobs = sum(is.na(a)),
                     ynobs_prop = ynobs/nd,
                     yobs_mean = mean(a,na.rm = T),
                     ex.under = case_when(
                       ynobs_prop==1 ~ y.min,
                       TRUE ~ y.min*ynobs_prop + yobs_mean *  (1-ynobs_prop)),
                     ex.upper = case_when(
                       ynobs_prop==1 ~ y.max,
                       TRUE ~ y.max*ynobs_prop + yobs_mean *  (1-ynobs_prop)) ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(prop.x = nd/sum(nd))
  
  nx.unique = dim(x.unique)[1] #number of unique x
  
  
  # * to compute: solve(\int xtx dFx) --------------------------------------------------------
  x.unique.list = split(x.unique[,1:nx],seq(nrow(x.unique[,1:nx])))
  xtx.unique = lapply(x.unique.list,function(h) (as.numeric(h) %*% t(as.numeric(h))) )
  xtx.final = Reduce('+',lapply(1:nx.unique, function(h) xtx.unique[[h]]*x.unique$prop.x[h]))
  
  xtx.m1 <- solve(xtx.final)
  
  
  # * to compute: gbar under and upper --------------------------------------
  x.unique$alpha.xt =  unlist(lapply(x.unique.list,function(h) as.numeric(cd%*% xtx.m1 %*% as.numeric(h))))
  
  x.unique = x.unique %>% 
    dplyr::mutate(
      gbar.under = case_when(
        alpha.xt > 0 ~ ex.under,
        TRUE ~ ex.upper),
      gbar.upper = case_when(
        alpha.xt > 0 ~ ex.upper,
        TRUE ~ ex.under),
      gbar.under.prop = gbar.under * prop.x,
      gbar.upper.prop = gbar.upper * prop.x
      )
  
  # * to compute: \int xt gbar_under/upper dFx ------------------------------
  xt.gunder <- Reduce('+',lapply(1:nx.unique,function(h) x.unique.list[[h]] * x.unique$gbar.under.prop[h]))
  xt.gupper <- Reduce('+',lapply(1:nx.unique,function(h) x.unique.list[[h]] * x.unique$gbar.upper.prop[h]))
  
  
  # * limits ----------------------------------------------------------------
  lb = cd%*% xtx.m1 %*%as.numeric(xt.gunder)
  ub = cd%*% xtx.m1 %*%as.numeric(xt.gupper)
  
  interval = c(lb,ub)
  names(interval) <- paste(x.names.or[-1][which(cd==1)],c("Li","Ls"))
  
  return(interval)
}

# APPLICATION -------------------------------------------------------------

# *data -------------------------------------------------------------------

chile.data <- rio::import("data/cs_biologicas.txt")
names(chile.data)

chile.data.x <- chile.data %>% 
  dplyr::select("PGA_1ERsemestre","matem","ptje nem","ranking","leng y Com")


# * descriptive statistics --------------------------------------------------
# total data
chile.data <- chile.data %>% 
  dplyr::mutate(selected = case_when(
    !is.na(PGA_1ERsemestre) ~ "Selected",
    TRUE ~ "Non selected"
  ))

n.stud <- dim(chile.data)[1]
n.stud

porc.by_sel_sexo <- chile.data %>% 
  dplyr::group_by(selected,sexo) %>% 
  dplyr::summarise(nd = n(),
                   mean.gpa = mean(PGA_1ERsemestre,na.rm=T),
                   sd.gpa = sd(PGA_1ERsemestre,na.rm=T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(selected) %>% 
  dplyr::mutate(n.sel = sum(nd)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(porc.of_sel = round(nd/n.sel,4)*100,
                porc.of_total = round(nd/sum(nd),4)*100,
                porc.bysel_oftotal = round(n.sel/sum(nd),4)*100) %>% 
  dplyr::ungroup() 
porc.by_sel_sexo

stat_factors_bysel <- chile.data %>% 
  dplyr::group_by(selected) %>% 
  dplyr::summarise(min_mat = min(matem,na.rm = T),
                   med_mat = median(matem,na.rm=T),
                   mea_mat = mean(matem,na.rm=T),
                   max_mat = max(matem,na.rm = T),
                   sdd_mat = sd(matem,na.rm=T),
                   min_leng = min(`leng y Com`,na.rm = T),
                   med_leng = median(`leng y Com`,na.rm=T),
                   mea_leng = mean(`leng y Com`,na.rm=T),
                   max_leng = max(`leng y Com`,na.rm = T),
                   sdd_leng = sd(`leng y Com`,na.rm=T),
                   min_rank = min(ranking,na.rm = T),
                   med_rank = median(ranking,na.rm=T),
                   mea_rank = mean(ranking,na.rm=T),
                   max_rank = max(ranking,na.rm = T),
                   sdd_rank = sd(ranking,na.rm=T),
                   min_nem = min(`ptje nem`,na.rm = T),
                   med_nem = median(`ptje nem`,na.rm=T),
                   mea_nem = mean(`ptje nem`,na.rm=T),
                   max_nem = max(`ptje nem`,na.rm = T),
                   sdd_nem = sd(`ptje nem`,na.rm=T)
                   )
stat_factors_table <- stat_factors_bysel %>% 
  tidyr::pivot_longer(cols=contains(c("_mat","_len","_rank","_nem")),names_to = "stat",values_to = "stat_val") %>% 
  dplyr::mutate(statt=substr(stat,start = 1,stop = 3),
                factor = substr(stat,start = 5,stop = nchar(stat))) %>% 
  dplyr::select(-stat) %>% 
  tidyr::pivot_wider(names_from = statt,values_from = stat_val)
  

stat_factors_table_sumary <- inner_join(stat_factors_table %>% dplyr::filter(selected=="Selected") %>% dplyr::select(-selected),
           stat_factors_table %>% dplyr::filter(selected=="Non selected") %>% dplyr::select(-selected),
           by="factor",suffix=c("_sel","_nonsel"))
  


# * plots -----------------------------------------------------------------

chile.data_longer = chile.data %>%
  dplyr::rename(leng_f=`leng y Com`,
                matem_f=matem,
                rank_f=ranking,
                nem_f=`ptje nem`) %>% 
  dplyr::select(id,selected,contains("_f")) %>% 
  tidyr::pivot_longer( cols = contains("_f"),values_to = "factor_val", names_to = "factor" ) 

chile.data_longer$factor <- factor(chile.data_longer$factor,
                                   levels = c("matem_f","leng_f","rank_f","nem_f"))
chile.data_longer <- chile.data_longer %>% 
  dplyr::arrange(factor)
  
boxplot_factors <- chile.data_longer %>% 
  ggplot(aes(x=factor,y=factor_val,fill=selected)) +
  geom_boxplot() +
  scale_fill_manual(values = c("violet","skyblue"),labels=c("Non enrolled","Enrolled")) +
  scale_x_discrete(labels=c("Mathematics","Language","Ranking","HS-GPA")) +
  guides(fill=guide_legend(title = "Applicant status"),) +
  xlab("Selection factor") +
  ylab("Score") + 
  theme_light() +
  theme(legend.position="bottom") +
  labs(title = "(a)")
  
selected <- chile.data %>% 
  dplyr::filter(selected=="Selected")


pga.summary <- summary(selected$PGA_1ERsemestre)
pga.mean <- mean(selected$PGA_1ERsemestre)
pga.sd <- sd(selected$PGA_1ERsemestre)
pga.min <- min(selected$PGA_1ERsemestre)
pga.max <- max(selected$PGA_1ERsemestre)

gpa_observed <- selected %>% 
  ggplot(aes(y=PGA_1ERsemestre,fill=selected)) +
  geom_boxplot() +
  scale_fill_manual(values = c("skyblue"))+
  ylab("GPA") + 
  theme_light() +
  theme(axis.text.x=element_blank(),legend.position = "none") +
  ylim(1,7)+
  labs(title = "(b)")
  
plot_sores_gpa <- grid.arrange(boxplot_factors, gpa_observed, ncol=2)

# * x complete -------------------------------------------------------
chile.data.xcomplete = chile.data.x[complete.cases(chile.data.x[,-1]),]

chile.data.xcomplete= chile.data.xcomplete %>% 
  dplyr::mutate( 
    PGA_1ERsemestre = (PGA_1ERsemestre - mean(PGA_1ERsemestre,na.rm=T))/sd(PGA_1ERsemestre,na.rm = T),
    matem = (matem - mean(matem))/sd(matem),
    `ptje nem` = (`ptje nem` -mean(`ptje nem` ))/sd(`ptje nem` ),
    ranking = (ranking - mean(ranking))/sd(ranking),
    `leng y Com` = (`leng y Com` - mean(`leng y Com`))/sd(`leng y Com`))
  

stoye_mat <- stoye_limits(datos = chile.data.xcomplete,Ymin=(1.0 - pga.mean)/pga.sd,Ymax = (7.0 - pga.mean)/pga.sd,cd=c(1,0,0,0,0))
stoye_nem <- stoye_limits(datos = chile.data.xcomplete,Ymin=(1.0 - pga.mean)/pga.sd,Ymax = (7.0 - pga.mean)/pga.sd,cd=c(0,1,0,0,0))
stoye_rank <- stoye_limits(datos = chile.data.xcomplete,Ymin=(1.0 - pga.mean)/pga.sd,Ymax = (7.0 - pga.mean)/pga.sd,cd=c(0,0,1,0,0))
stoye_leng <-  stoye_limits(datos = chile.data.xcomplete,Ymin=(1.0 - pga.mean)/pga.sd,Ymax = (7.0 - pga.mean)/pga.sd,cd=c(0,0,0,1,0))
stoye_inter <- stoye_limits(datos = chile.data.xcomplete,Ymin=(1.0 - pga.mean)/pga.sd,Ymax = (7.0 - pga.mean)/pga.sd,cd=c(0,0,0,0,1))


stoye_table <- rbind(
      stoye_mat,
      stoye_leng,
      stoye_rank,
      stoye_nem)
colnames(stoye_table) <- c("LB","UB")
rownames(stoye_table) <- NULL


stoye_table <- stoye_table %>% 
  as_tibble() %>% 
  dplyr::mutate(factor=c("Math","Language","Ranking","HS-GPA")) %>% 
  dplyr::select(factor,LB,UB) %>% 
  dplyr::mutate(across(where(is.numeric),round,4))

stoye_mat_minmax <- stoye_limits(datos = chile.data.xcomplete,Ymin=(pga.min - pga.mean)/pga.sd,Ymax = (pga.max - pga.mean)/pga.sd,cd=c(1,0,0,0,0))
stoye_nem_minmax <- stoye_limits(datos = chile.data.xcomplete,Ymin=(pga.min - pga.mean)/pga.sd,Ymax = (pga.max - pga.mean)/pga.sd,cd=c(0,1,0,0,0))
stoye_rank_minmax <- stoye_limits(datos = chile.data.xcomplete,Ymin=(pga.min - pga.mean)/pga.sd,Ymax = (pga.max - pga.mean)/pga.sd,cd=c(0,0,1,0,0))
stoye_leng_minmax <-  stoye_limits(datos = chile.data.xcomplete,Ymin=(pga.min - pga.mean)/pga.sd,Ymax = (pga.max - pga.mean)/pga.sd,cd=c(0,0,0,1,0))
stoye_inter_minmax <- stoye_limits(datos = chile.data.xcomplete,Ymin=(pga.min - pga.mean)/pga.sd,Ymax = (pga.max - pga.mean)/pga.sd,cd=c(0,0,0,0,1))


stoye_table_minmax <- rbind(
  stoye_mat_minmax,
  stoye_leng_minmax,
  stoye_rank_minmax,
  stoye_nem_minmax)
colnames(stoye_table_minmax) <- c("LB","UB")
rownames(stoye_table_minmax) <- NULL


stoye_table_minmax <- stoye_table_minmax %>% 
  as_tibble() %>% 
  dplyr::mutate(factor=c("Math","Language","Ranking","HS-GPA")) %>% 
  dplyr::select(factor,LB,UB) %>% 
  dplyr::mutate(across(where(is.numeric),round,4))

# * xy complete -----------------------------------------------------------
chile.data.xycomplete = chile.data.x[complete.cases(chile.data.x),]

chile.data.xycomplete = chile.data.xycomplete %>% 
  dplyr::mutate( 
  PGA_1ERsemestre = (PGA_1ERsemestre - mean(PGA_1ERsemestre,na.rm=T))/sd(PGA_1ERsemestre,na.rm = T),
  matem = (matem - mean(matem))/sd(matem),
  `ptje nem` = (`ptje nem` -mean(`ptje nem` ))/sd(`ptje nem` ),
  ranking = (ranking - mean(ranking))/sd(ranking),
  `leng y Com` = (`leng y Com` - mean(`leng y Com`))/sd(`leng y Com`))

anova_factors = lm(chile.data.xycomplete$PGA_1ERsemestre ~ chile.data.xycomplete$matem + chile.data.xycomplete$`ptje nem` + chile.data.xycomplete$ranking +chile.data.xycomplete$`leng y Com`)
anova_factors
anova(anova_factors)
anova_factors_table <- tibble(factor=c("Math","NEM","Ranking","Language"),
      value= round(as.numeric(coefficients(anova_factors)[-1]),4))

anova_factors_table$factor <- factor(anova_factors_table$factor,
                                     levels = c("Math","Language","Ranking","NEM"))

anova_factors_table = anova_factors_table %>% 
  dplyr::arrange(factor)
