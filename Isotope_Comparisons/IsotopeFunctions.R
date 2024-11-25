


### CV of different fish groups 
# ---------------------- Library Loads ---------------------------------------
library(RColorBrewer)
library(abind)

#----------------------------- Params ---------------------------------------- 

x.p = seq(-38,-20,length.out=100); 
y.p = seq(0,9, length.out = 100) 
parms <- list();
parms$n.iter=2*10^4; 
parms$n.burnin=1*10^3; 
parms$n.thin=10
parms$n.chains=3     
priors=list();
priors$R=1*diag(2);
priors$k=2; 
priors$tau.mu=1.0E-3 # Vague priors 
Nsamples=1000
n.posts <- 1000;
p.ell <- 0.90 # How much data to include? Standard ellipses --> p.ell = .9
n.points = 1000




overlap = function(data_input, comm, dr, posterior){
  ## Remove this later
  data_overlap = data_input %>% 
    filter(community == comm)
  
  #data_overlap = combo %>% filter(community == 1)
  spp=length(unique(data_overlap$group))
  print(data_overlap)
  
  name_list = names(posterior) 
  
  
  name_matrix = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","G1"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F) %>%
    filter(C1 == comm & C2 == comm) %>%
    select(Var1, Var2) %>%
    filter(Var1 != Var2)
  
  
  names_modified = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","group"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F)%>%
    filter(C1 == comm & C2 == comm) %>%
    mutate(group = as.numeric(group))%>%
    left_join(legend %>% 
                mutate(community = as.character(community)), by = c("group", "C1" = "community")) %>%
    select(-color, -common) %>%
    rename("G1" = group) %>%
    rename(group = "G2") %>%
    rename(Code1 = CODE) %>%
    mutate(group = as.numeric(group)) %>%
    left_join(legend %>% 
                mutate(community = as.character(community)), by = c("group", "C1" = "community")) %>%
    select(-color, -common) %>%
    rename(G2 = group) %>%
    unite(col = Name, c(Code1, CODE), sep = " v " ) %>%
    filter(Var1 != Var2) %>%
    select(Name)
  
  
  ma = matrix(NA, nrow = dr, ncol = length(name_matrix[,1]))
  
  for(i in 1:length(name_matrix[,1])){ 
    
    overlap = bayesianOverlap(ellipse1 = as.character(name_matrix[i,1]), 
                              ellipse2 = as.character(name_matrix[i,2]), 
                              posterior, 
                              draws = dr,
                              p.interval = 0.95,
                              n = dr)
    
    bayes.prop <- (overlap[,3] / (overlap[,2] + overlap[,1] - overlap[,3]))
    
    ma[,i] = bayes.prop
    
  }
  
  
  
  full_olap_mat = ma %>% as.data.frame() %>%
    rename_with(~ setNames((unique(names_modified[,1])), .), everything())
  
  full_olap_mat[dr+1,] = unique(names_modified[,1])
  
  rownames(full_olap_mat) = c((data.frame(n = "com",c = comm, post = c(1:dr)) %>% unite("rowname", c(n, c, post)))$rowname, "SppPair")
  
  
  #return(as.data.frame(t(olap_mat))) ## This is a reminder to pull 90% CI intervals out of this at a point when i can play with code and don't need to be writing
  return(t(as.data.frame(full_olap_mat)))
  
}


ellip_data = function(numb_species, numb_posts, posterior){
  for(i in 1:numb_species){
    dat = vector()
    for(j in 1:numb_posts){
      ellipse_data =  ellipse::ellipse(x = matrix(c(posterior[[2]][[i]][j,1],
                                                    median(posterior[[2]][[i]][j,2]),
                                                    median(posterior[[2]][[i]][j,3]),
                                                    median(posterior[[2]][[i]][j,4])),
                                                  2,2), 
                                       centre = c(median(posterior[[2]][[i]][j,5]),
                                                  median(posterior[[2]][[i]][j,6])), level = .95, 
                                       npoints = n.points) %>% as.data.frame() %>%
        summarise_all(list("min"=min, "max"=max)) %>% as.matrix()
      
      dat = rbind(dat, ellipse_data)
    }
    
    ellip[,,i] = dat
    
  }
  
  return(ellip)
}
