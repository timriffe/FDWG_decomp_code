# remotes::install_github("timriffe/LifeIneqDecomp")
# remotes::install_github("timriffe/DemoDecomp",force=TRUE)
library(HMDHFDplus)
library(LifeIneqDecomp)
library(LifeIneq)
library(tidyverse)
library(DemoDecomp)



Mx <- readHMDweb("ESP","Mx_1x1",username = us, password = pw)%>% 
  filter(Year %in% c(1950,2000)) %>% 
  select(-OpenInterval) %>% 
  pivot_longer(Female:Total, names_to = "Sex", values_to = "Mx") %>% 
  mutate(Mx = ifelse(is.na(Mx),0,Mx))
Ex <- readHMDweb("ESP","Exposures_1x1",username = us, password = pw)%>% 
  filter(Year %in% c(1950,2000))%>% 
  select(-OpenInterval) %>% 
  pivot_longer(Female:Total, names_to = "Sex", values_to = "Px")


Dat <- Mx %>% 
  left_join(Ex, by = c("Year","Age","Sex")) %>% 
  pivot_wider(names_from = Year, values_from = c(Mx, Px)) %>% 
  arrange(Sex, Age) %>% 
  select(Sex, Age, 3:6)

write_csv(Dat, file = "example_data.csv")



Mx %>% 
  filter(Year %in% c(1950, 2018)) %>% 
  ggplot(aes(x = Age, y = Total, color = as.factor(Year), group = Year)) +
  geom_line() + 
  scale_y_log10()

# some low-effort lifetable functions:
mx2lx <- function(mx){
  mx[is.na(mx)] <- 0
  lx <- exp(-cumsum(mx))
  lx <- c(1,lx)
  lx[1:length(mx)]
}

lx2dx <- function(lx){
  -diff(c(lx,0))
}
  
lx2ex <- function(lx){
  Lx <- (lx + c(lx[-1],0)) / 2
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex
}

the_goods <- 
Mx %>% 
  filter(Year %in% c(1950,2018)) %>% 
  select(Year, Age, mx = Total) %>% 
  group_by(Year) %>% 
  mutate(lx = mx2lx(mx),
         dx = lx2dx(lx),
         ex = lx2ex(lx),
         ax = .5) %>% 
  ungroup() %>%
  select(-mx) %>% 
  pivot_wider(names_from = Year, values_from = c(ax,dx,lx,ex))

ax <- the_goods %>% select(contains("ax")) %>% as.matrix()
lx <- the_goods %>% select(contains("lx")) %>% as.matrix()
dx <- the_goods %>% select(contains("dx")) %>% as.matrix()
ex <- the_goods %>% select(contains("ex")) %>% as.matrix()

bw_decomp(age = 0:110, 
          ax,
          dx,
          lx,
          ex,
          prop = c(.5, .5),
          method = "gini")

# test

my_gini <- function(mx){
  
    age      <- seq_along(mx) - 1
    N        <- length(mx)
    # derive what we need
    ax       <- rep(.5,N)
    lx       <- mx2lx(mx)
    dx       <- lx2dx(lx)
    ex       <- lx2ex(lx)
    
    # matrix of age differences
    ad       <- outer(age + ax, age + ax, "-")
    
    # joint probability of death
    pd       <- outer(dx, dx, "*")

    sum(abs(ad) * pd) / (2 * ex[1])
}


# mx2gini <- function(mx, age_interval = 1){
#   age <- (seq_along(mx) - 1) * age_interval
#   lx <- mx2lx(mx)
#   ex <- lx2ex(lx)
#   ax <- rep(.5, length(mx))
#   ineq_gini(age=age, ax = ax, lx = lx, ex = ex)
# }

# Which ages contribute to difference in gini?
library(magrittr)
Mx %>% 
  filter(Year %in% c(1950,2000)) %>% 
  select(Age, mx = Total) 
  

mx <-
  Mx %>% 
  dplyr::filter(Year == 2000) %>% 
  dplyr::pull(Total)

age <- 0:110

my_gini(mx)

Mx %>% 
  filter(Year %in% c(2000,2010)) %>% 
  select(Year,Age, mx = Total) %>% 
  pivot_wider(names_from = Year, values_from = mx) %>% 
  mutate(ho = horiuchi(my_gini, `2000`, `2010`, N = 20),
         ca = ltre(my_gini, `2000`, `2010`),
         sp = stepwise_replacement(my_gini, `2000`, `2010`)) %>% 
  select(Age, ho, ca, sp) %>% 
  pivot_longer(2:4, names_to = "variant", values_to = "cx") %>% 
  ggplot(aes(x = Age, y = cx, color = variant, lty = variant)) + 
  geom_line()

# plot(numDeriv::grad(my_gini,mx))
# plot(gini_dot(mx1,mx2))
# plot(gini_dot(mx1,mx2) * (mx1 - mx2))
# 
# gini_dot <- function(mx1, mx2){
#  mx <- (mx1 + mx2)/2
#  rhox <- (mx1 - mx2) / mx
#  lx <- mx2lx(mx)
#  ex <- lx2ex(mx)
#  
#  D  <- 1 - my_gini(mx)
#  Dx <- rev(cumsum(rev(lx ^2))) / rev(cumsum(rev(lx))) * (1 / lx)
#  wx <- mx * lx * ex
#  Wx <- (2 * lx * Dx - D) / ex[1]
#  
#  rhox * wx * Wx
# }
# 
# mx1 <- Mx %>% 
#   filter(Year == 1950) %>% pull(Total)
# mx2<- Mx %>% 
#   filter(Year == 2000) %>% pull(Total)

# cumsum(mx)
# age[cumsum(mx) >= 1]
# mx2lx(mx) %>% lx2ex()



