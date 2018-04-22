# install.packages("igraph")# if this failed
# install.packages("igraph", type="binary")
# install.packages("pracma", type="binary")

library(igraph)
library(Matrix)
library(pracma)

set.seed(1)

create_transition_matrix = function (g){
  
  # WARNING: make sure your graph is connected (you might input GCC of your graph)
  
  vs = V(g)
  n = vcount(g)
  adj = as_adjacency_matrix(g)
  adj[diag(rowSums(adj) == 0)] = 1  # handle if the user is using the function for networks with isolated nodes by creating self-edges
  z = matrix(rowSums(adj, , 1))
  
  transition_matrix = adj / repmat(z, 1, n)  # normalize to get probabilities
  
  return(transition_matrix)
}

# random_walk = function (g, num_steps, start_node, transition_matrix = NULL){
#   if(is.null(transition_matrix)){
#     transition_matrix = create_transition_matrix(g)
#   }
#   
#   v = start_node
#   for(i in 1:num_steps){
#     #         fprintf('Step %d: %d\n', i, v)  # COMMENT THIS
#     PMF = transition_matrix[v, ]
#     v = sample(1:vcount(g), 1, prob = PMF)        
#   }
#   
#   return(v)
# }

random_walk_plot = function(rn) {
  avg_dis_all <- numeric()
  avg_std_all <- numeric()
  
  degreesVector <- degree(rn)
  degreesEnd <- numeric()
  
  num_node <- vcount(rn)
  
  
  for(step in 1:20){
    dis <- numeric()
    avg_dis <- numeric()
    avg_std <- numeric()
    
    if(num_node > 1000){
      sample = sample(vcount(rn),1000)
      sample
      for (start_node in sample) {
        v_last = tail(random_walk(rn, start_node, step), n=1)
        dis <- c(dis, shortest.paths(rn, start_node, v_last))
        
        degreesEnd <- c(degreesEnd, degreesVector[v_last])
      }
    }
    else{
      for (start_node in 1:num_node) {
        v_last = tail(random_walk(rn, start_node, step), n=1)
        dis <- c(dis, shortest.paths(rn, start_node, v_last))
        
        degreesEnd <- c(degreesEnd, degreesVector[v_last])
      }
    }
    
    
    avg_dis <- mean(dis)
    avg_std <- mean((dis - mean(dis))**2)
    
    fprintf('Step %d: %f %f\n', step, avg_dis, avg_std)  
    avg_dis_all <- c(avg_dis_all, avg_dis)
    avg_std_all <- c(avg_std_all, avg_std)
    
  }
  
  #plots
  title = "Erdös-Rényi"
  
  filename = paste('fig/','avgdis_', title, '_', num_node,'.pdf', sep="")
  plot(avg_dis_all, typ='l', main = paste("Average Distance v.s. t,", title, num_node, "nodes"), xlab = "Steps", ylab = "Average Distance")

  
  plot(avg_std_all, typ='l', main = paste("Average Standard Deviation v.s. t,", title, num_node, "nodes"), xlab = "Steps", ylab = "Average Standard Deviation")
  
  degreesVector <- degree(rn)
  
 
  hist(degreesVector, main = "Degree distribution of graph", xlab="Degrees")
  
  
  hist(degreesEnd, main = "Degree distribution of the nodes reached at the end", xlab="Degrees")
}

rn <- erdos.renyi.game(1000, 0.01, type="gnp")
random_walk_plot(rn)

rn <- erdos.renyi.game(100, 0.01, type="gnp")
random_walk_plot(rn)

rn <- erdos.renyi.game(10000, 0.01, type="gnp")
random_walk_plot(rn)