# Use 15N data to estimate the range of cannibalism in the food web presented by Oelbermann and Scheu 2010.

# Author: Robert W. Buchkowski
# Date Created: Oct. 6, 2021

# Load in libraries:
library(pacman)
p_load(tidyverse)

# Oelbermann and Scheu 2010 ----

# Scrape data from the paper:
if(F){ # Only run once. Data saved in the directory
  p_load(tabulizer)
  
  Oel10 <- extract_areas("C:/Users/rober/Zotero/storage/TMYI46QK/Oelbermann and Scheu - 2010 - Trophic guilds of generalist feeders in soil anima.pdf", pages = c(9,10))
  
  Oel10 <- do.call("rbind", Oel10)
  
  Oel10 <- as.data.frame(Oel10)
  colnames(Oel10) <- c("TL", "Pool", "x15N", "Habitat")
  
  Oel10 = Oel10 %>% tibble() %>%
    mutate(x15N = str_replace(x15N, "x","-")) %>%
    separate(x15N, into = c("x15N", "sd15N"), sep = "[+]") %>%
    mutate(x15N = as.numeric(x15N),
           sd15N = as.numeric(sd15N))
  
  Oel10 %>% write_rds("Data/Oelbermann_Scheu_2010.rds")
}

# Run the analysis
Oel10 <- read_rds("Data/Oelbermann_Scheu_2010.rds")

Oel10 %>%
  filter(TL != "") %>%
  select(-Pool, -Habitat) %>%
  ggplot(aes(x = TL, y = x15N, ymin = x15N-sd15N, ymax = x15N+sd15N)) + geom_pointrange()

# next task is to assign feeding relationships and test the results by trophic group.

Oelmat <- matrix(0, nrow = 18, ncol = 18, dimnames = list(unique(Oel10$TL)[-2],unique(Oel10$TL)[-2]))

# Row eats column
Oelmat[c("D1", "D2", "D3", "D4", "D5", "D6", "D7"), "B"] = 1

Oelmat[c("P1", "P2", "P3", "P4"), c("H1", "H2", "H3", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "P1", "P2", "P3", "P4")] = 1

Oelmat[c("N1", "N2", "N3"), c("H1", "H2", "H3", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "P1", "P2", "P3", "P4","N1", "N2", "N3")] = 1

# Convert feeding matrix to Predator-Prey list:
Oelmat <- as.data.frame(Oelmat)

Oelmat[,"Predator"] = row.names(Oelmat)
row.names(Oelmat) = NULL

Oelmat = Oelmat %>% pivot_longer(-Predator,names_to = "Prey")

Oel15N = Oel10 %>%
  filter(TL != "") %>%
  select(-Pool, -Habitat)

Oelmat %>%
  filter(value == 1) %>%
  select(-value) %>%
  left_join(
    Oel15N %>% rename(Predx15N = x15N, Predsd15N = sd15N), by = c("Predator" = "TL")
  ) %>%
  left_join(
    Oel15N, by = c("Prey" = "TL")
  ) %>%
  # Calculation with no uncertainty:
  mutate(d15N = Predx15N - x15N) %>%
  # remove canniablism:
  filter(Predator != Prey) %>%
  ggplot(aes(x = Predator, y = d15N, color = Prey)) + geom_hline(yintercept = 3.4) + stat_summary(aes(group = Predator),fun = mean, color = "black") + geom_point()

# Run a constraint model like the one in the Frontiers paper to see how much a predator can eat from cannibalism and still get the 15N value seen in the paper.

# Try to calculate the mixing model:

for(i in 8:11){
  f.obj = rep(0,11)
  f.obj[i] = 1
  f.con = matrix(0,ncol = 11, nrow = 11)
  diag(f.con) = 1
  f.con = rbind(Oel15N$x15N[5:15]+3.4,rep(1, 11),f.con)
  f.dir = c("=", "=", rep(">=", 11))
  f.rhs = c(Oel15N$x15N[5:15][i],1, rep(0,11))
  
  outputU = lpSolve::lp("max",
                        f.obj,
                        f.con,
                        f.dir,
                        f.rhs)
  
  outputL = lpSolve::lp("min",
                        f.obj,
                        f.con,
                        f.dir,
                        f.rhs)
  
  print(c(outputL$solution[i], outputU$solution[i]))
}



# Jassey et al. 2013----

# Scrape data from the paper:
if(F){ # Only run once. Data saved in the directory
  p_load(tabulizer)
  
  J13 <- extract_areas("C:/Users/rober/Zotero/storage/IDG8IH4U/Jassey et al. - 2013 - To What Extent Do Food Preferences Explain the Tro.pdf", pages = c(4))
  
  J13 <- J13[[1]]
  
  write_rds(J13, "Data/J13.rds")
  
  J13a = J13[3:18,c(2,5,6)]
  
  J13a = J13a[J13a[,2] != "",]
  J13a[c(2,3,10,11),1] = c("POMB", "Microalgae", "Rotifers", "Nematodes")
  colnames(J13a) = c("ID", "13C", "15N")
  
  J13a = tibble(data.frame(J13a)) %>%
    separate(X13C, into = c("X13C", NA), sep = "±") %>%
    separate(X15N, into = c("X15N", NA), sep = "±") %>%
    separate(X13C, into = c(NA,"X13C"), sep = 1) %>%
    separate(X15N, into = c(NA,"X15N"), sep = 1) 
  
  J13a[c(5,11),3] = c("1.0", "0.5")
  
  J13a = J13a %>%
    mutate(X13C = as.numeric(X13C)*-1) %>%
    mutate(X15N  = as.numeric(X15N)*c(-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1))
  
  write_rds(J13a, "Data/Jassey2013_isotope.rds")
}

# Load in the data
J13a = read_rds("Data/Jassey2013_isotope.rds")

# Max or min cannibalism
f.obj = c(rep(0,7),1,0,0,0)

# Set the constraint matrix
f.con = matrix(0,ncol = 11, nrow = 11)
diag(f.con) = 1
f.con = rbind(J13a$X15N+3.4,rep(1, 11),f.con)

# Set the direction
f.dir = c("=", "=", rep(">=", 11))

# Set the rhs
f.rhs = c(J13a$X15N[8],1, rep(0,11))

outputU = lpSolve::lp("max",
                      f.obj,
                      f.con,
                      f.dir,
                      f.rhs)

outputL = lpSolve::lp("min",
                      f.obj,
                      f.con,
                      f.dir,
                      f.rhs)

# Range of possible rates:
c(outputL$solution[8],outputU$solution[8])
