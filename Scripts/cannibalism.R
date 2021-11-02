# Cannibalism analysis:
# Author: Robert Buchkowski
# Date created: July 17/2021

# Install soilfoodwebs from Github:

if(!require("soilfoodwebs")){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse,cowplot,gridGraphics, soilfoodwebs)

# Maximum cannibalism plot -----
# Try to plot the function for maximum cannibalism from the MS:
expand.grid(CEp = seq(0,0.5,by=0.01),
            CEc = seq(0,0.5, by = 0.01)) %>%
  tibble() %>%
  mutate(mcr = CEp/(1-CEc+CEp)) %>%
  ggplot(aes(x = CEc, y = CEp, fill = mcr)) + 
  geom_raster() + theme_classic() + 
  scale_fill_viridis_c(name = "Max cannibalism rate") +
  xlab("Conversion efficiency canibalism") +
  ylab("Conversion efficiency other prey") + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "top")


# Cannibalism ----
# Equation 5:
canrateopt <- function(canrate, COMM, node = 1){
  fedvec = COMM$imat[node,]
  fedvec[node] = canrate
  
  return((fedvec[node]*COMM$prop$B[node]/sum(fedvec*COMM$prop$B) - COMM$prop$a[node]*COMM$prop$p[node])^2)
}

cantable <- function(COMMUN){
  results <- rep(0, length = dim(COMMUN$imat)[1])
  
  for(selnode in 1:length(results)){
    res <- optimise(canrateopt, lower = 0, upper = 1000000000, COMM = COMMUN, node = selnode)
    if(res$objective < 1e-6){
      results[selnode] = res$minimum
    }else{
      results[selnode] = NA
      print(selnode)
    }
    
  }
  
  results[COMMUN$prop$isDetritus == 1 | is.na(results) | TLcheddar(COMMUN$imat) == 1] = 0
  
  imat_mod = COMMUN$imat
  diag(imat_mod) = results
  
  COMMUN2 = COMMUN
  
  COMMUN2$imat = imat_mod
  
  ji = Jacobsindex(COMMUN2)
  
  tibble(ID = COMMUN$prop$ID,
         MaxCan = COMMUN$prop$a*COMMUN$prop$p,
         PrefCan = signif(results/rowSums(imat_mod),3)) %>%
    full_join(
      tibble(ji) %>%
        filter(Predator == Prey) %>%
        select(-Prey) %>% rename(ID = Predator), by = "ID"
    )
}

# Analysis of all communities ----

# Create a table of the trophic species for each web and name groups:

if(F) { # Not used after initial creation, file modified in excel
  data.frame(ID = c(Koltz2018$prop$ID,
                    Holtkamp2011$Young$prop$ID,
                    Andres2016$GA$prop$ID,
                    CPER$prop$ID,
                    deRuiter1994$INT$prop$ID)) %>%
    write_csv("Data/Key.csv")
}

nametocolumn <- function(LIT){
  for(i in 1:length(LIT)){
    LIT[[i]] = cbind(LIT[[i]], Web = names(LIT)[i])
  }
  return(LIT)
}

result = rbind(
  do.call("rbind",nametocolumn(lapply(Andres2016, cantable))),
  do.call("rbind",nametocolumn(lapply(Holtkamp2011, cantable))),
  do.call("rbind",nametocolumn(lapply(deRuiter1994, cantable))),
  cbind(cantable(CPER), Web = "CPER"),
  cbind(cantable(Koltz2018), Web = "Koltz2018")
)

tibble(result) %>%
  filter(Jacobs < 0)

result %>%
  left_join(
    read_csv("Data/Key.csv"), by = "ID"
  ) %>%
  filter(!is.na(Group)) %>%
  pivot_longer(MaxCan:Jacobs) %>%
  distinct() %>%
  group_by(Group, name) %>%
  summarize(L = min(value, na.rm = T), U = max(value, na.rm = T), N = n()) %>%
  mutate(Res = paste0(signif(L,2)," to ",signif(U,2))) %>%
  select(Group, name, N, Res) %>%
  pivot_wider(values_from = Res)  %>% 
  write_excel_csv("Data/Cannibalism/table1.csv")

# Change in web efficiency ----
# This analysis will produce estimates of the change in web efficiency with and without cannibalism.

getcanndiff <- function(cann_rate = 0.99, testcomm, cancan, rtnall = F){
  caninput = cantable(testcomm) %>%
    select(ID, MaxCan)
  
  diag(testcomm$imat) = 0
  
  if(any(colnames(testcomm$imat) != caninput$ID)) stop("Erro matching tables.")
  
  caninput = caninput$MaxCan*cann_rate
  
  canrates = rowSums(testcomm$imat * matrix(testcomm$prop$B, ncol = dim(testcomm$imat)[1], nrow = dim(testcomm$imat)[1], byrow = T)) * caninput /
    ((1-caninput)*testcomm$prop$B)
  
  canrates[is.na(canrates)] = 0
  canrates[cancan == 0] = 0
  
  diag(testcomm$imat) = canrates
  
  check = comana(testcomm, shuffleTL = F, eqmtolerance = 1e-5)$fmat
  
  # Check the rates to make sure they came out OK
  if(any(diag(check)/rowSums(check) - caninput > 0.001, na.rm = T)) {warning("Calculations of rate not coming out right.")}
  
  output = comana(testcomm, shuffleTL = F, eqmtolerance = 1e-5)
  
  if(any(output$Nmin[output$usin$prop$canIMM == 0] < 0)){
    warning("Cannibalism causing negative Nmin.")
  }
  
  if(rtnall){
    return(output)
  }else{
    if(all(
      output$consumption[TLcheddar(testcomm$imat) != 1] >= 0
      )
    ){
      return(c(cann_rate = cann_rate,
               Cprop = sum(output$Cmin)/sum(output$fmat),
               Nprop = sum(output$Nmin)/sum(output$Nfmat)))
    }else{
      return(c(cann_rate = cann_rate,
               Cprop = NA,
               Nprop = NA))
    }
  }
}

# Create a vector for the results:
results = vector(mode = "list", length = 18)

# CPER
results[[1]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = CPER, cancan = c(1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0)))), Web = "CPER")

# Holtkamp
results[[2]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = corrstoich(Holtkamp2011$Young), cancan = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)))), Web = "Holtkamp Young")

results[[3]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = corrstoich(Holtkamp2011$Mid), cancan = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)))), Web = "Holtkamp Mid")

results[[4]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = corrstoich(Holtkamp2011$Old), cancan = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)))), Web = "Holtkamp Old")

results[[5]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = corrstoich(Holtkamp2011$Heathland), cancan = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)))), Web = "Holtkamp Heathland")

# Andres
results[[6]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Andres2016$GA, cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 GA")

results[[7]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Andres2016$UGA, cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 UGA")

results[[8]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(Andres2016$UGB), cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 UGB")

results[[9]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Andres2016$GB, cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 GB")

results[[10]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(Andres2016$UGC), cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 UGC")

results[[11]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Andres2016$GC, cancan = c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1)))), Web = "Andres2016 GC")

# deRuiter
  results[[12]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(deRuiter1994$CON), cancan = c(1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0)))), Web = "deRuiter1994 CON")

results[[13]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(deRuiter1994$CON10), cancan = c(1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0)))), Web = "deRuiter1994 CON10")

results[[14]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(deRuiter1994$INT), cancan = c(1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0)))), Web = "deRuiter1994 INT")

results[[15]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = checkcomm(deRuiter1994$INT10), cancan = c(1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0)))), Web = "deRuiter1994 INT10")

# Koltz

# ... Prepare community
Koltztemp = checkcomm(Koltz2018)
diag(Koltztemp$imat) = 0 # Remove cannibalism

Koltztemp$prop[Koltztemp$prop$ID == "Diatoms",c("a", "p")] = 1 # Make diatoms efficient to get rid of their carbon mineralization

results[[16]] <- cbind(
  data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Koltztemp, 
                      cancan = c(0,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 0,1,1,0,0,
                                 0,0,0,0,0,1)
                        ))), Web = "Koltz Original")

# Add in no mutual feeders 
tocombine = can_mutfeed(Koltztemp)$prop %>%
  select(ID, MutualPred) %>%
  filter(!is.na(MutualPred))

Koltztemp2 = comtrosp(Koltztemp, selected = tocombine$ID, deleteCOMBOcannibal = T, allFEEDING1 = 1)

tocombine = can_mutfeed(Koltztemp2)$prop %>%
  select(ID, MutualPred) %>%
  filter(!is.na(MutualPred))

Koltztemp2 = comtrosp(Koltztemp2, selected = tocombine$ID, deleteCOMBOcannibal = T, allFEEDING1 = 1)

results[[17]] <- cbind(
  data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = Koltztemp2, 
                      cancan = c(0,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 1,1,1,1,1,
                                 0,1,1,0,0,
                                 0,0,0,0,0,0)
  ))), Web = "Koltz No mutual predation")

# Koltz in CPER web

# ... Create new community
if(F){ # Not used after the initial creation: Resulting file is saved in the directory
  
  KoltzCPER = CPER
  
  write_csv(CPER$prop, "Data/CPERprop2.csv")
  write_csv(Koltz2018$prop, "Data/Koltzprop.csv")
  
  newprop = read_csv("Data/CPERprop.csv")
  
  KoltzCPER$prop = KoltzCPER$prop %>%
    select(ID) %>%
    left_join(
      newprop, by = "ID"
    )

  KoltzCPER = checkcomm(KoltzCPER)
  KoltzCPER %>% write_rds("Data/KoltzCPER.rds")
}

KoltzCPER <- read_rds("Data/KoltzCPER.rds")
results[[18]] <- cbind(data.frame(t(sapply(seq(0, 0.99, by = 0.01), getcanndiff, testcomm = KoltzCPER, cancan = c(1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0)))), Web = "Koltz data in CPER")

# Final plot
resultsf = tibble(do.call("rbind", results))


resultsf %>%
  pivot_longer(Cprop:Nprop) %>%
  ggplot(aes(x = cann_rate, y = value, color = Web)) + geom_line() + facet_wrap(.~name, scales = "free") + theme_classic()


png("Figure2.png", width = 9, height = 4, units = "in", res = 600)
resultsf %>%
  mutate(Web2 = Web) %>%
  separate(Web2, into = c("Manuscript", NA, NA,NA,NA,NA), sep = " ") %>%
  rename(Carbon = Cprop, Nitrogen = Nprop) %>%
  pivot_longer(Carbon:Nitrogen) %>%
  left_join(
    tibble(Web = unique(resultsf$Web),
           Treatment = c("A", "A", "B", "C", "D",
                         "A", "B", "C", "D", "E", "F",
                         "A", "B", "C", "D",
                         "A", "Koltz: No mutual pred.", "Koltz in CPER"),
           Group     = factor(c("None", "None", "None", "None", "None",
                         "None", "None", "None", "None", "None",
                         "None", "None", "None", "None", "None",
                         "None", "Koltz: No mutual pred.", "Koltz in CPER"), levels = c("None", "Koltz: No mutual pred.", "Koltz in CPER"))), by = "Web"
  ) %>%
  ggplot(aes(x = cann_rate, y = value, color = Manuscript, linetype = Group, group = Web)) + geom_line() + facet_wrap(.~name, scales = "free") + theme_classic() + ylab("Mineralization (prop. of flux)") + xlab("Cannibalism rate (prop. of maximum)")
dev.off()