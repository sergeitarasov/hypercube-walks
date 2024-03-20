

Q <- initQ(c(0,1), c(1), diag.as = 0)

Qj <-  amaSMM(Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, controlling.state = NULL)

#Qj <-  amaSMM(Q, Q, Q,  controlling.state = NULL)
nrow(Qj)
# Convert Q matrix to adjacency matrix
#adj_matrix <- ifelse(Qj == 1, 1, 0)



# Create graph object
G <- graph.adjacency(Qj, mode = "undirected")

# Plot the graph
plot(G, 
     vertex.size = 30, 
     vertex.color = "lightblue", 
     vertex.label.cex = 1.5,
     edge.width = 2,
     edge.color = "gray")

# we leave the first state, and sample states to remove
s2remove=sample_states(total=1024, 400)
#s2remove=c(7,8)
Qr <- Qj[-s2remove,-s2remove]

#----


#get_neighbours(Qj, focal.state='111', max.dist=1)
s2remove= get_neighbours_maxdist(Qj, focal.state='1111111111', max.dist=5)
s2remove
length(s2remove)
Qr <- Qj[-s2remove,-s2remove]


s2remove2= get_neighbours_maxdist(Qr, focal.state='1111011111', max.dist=4)
s2remove2
length(s2remove2)
Qr <- Qr[-s2remove2,-s2remove2]

#--- mutate
initial_string='0000000000'
selected_states <- c()
for (i in 1:1){
  x <- mutate_binary_string(initial_string, generations=100)
  selected_states <- c(x, selected_states)
}
selected_states <- unique(selected_states)
selected_states
idx <- match(selected_states, rownames(Qj))
Qr <- Qj[idx,idx]

#----
Gr <- graph.adjacency(Qr, mode = "undirected")
plot(Gr, 
     vertex.size = 5, 
     edge.width = 2,
     vertex.label = NA,
     edge.color = "gray")
# Calculate connected components
comp <- components(Gr)
comp$no


Qscl <- adj2Q(Qr) 
get_lumpQrates_all(Qscl, nchar=10)



#------------ get approximate rate estimates
colnames(Qr)
Qscl
get_lumpQrates(Qscl,1)

Q.reod <- reoder_Q(Qscl, 2)
get_lumpQrates(Q.reod,2)

Q.reod <-reoder_Q(Qscl, 3)
get_lumpQrates(Q.reod,3)

get_lumpQrates_all(Qscl, nchar=3)

get_lumpQrates_all(Qscl, nchar=10)

#-- simulate
tree<-pbtree(n=100, scale=100, b=1, d=0)
plot(tree)
rowSums(Qscl)
Qscl[1,]
Qsim=Qscl*.1
hist <- sim.history(tree, Qsim, nsim=1, anc='0000000000')
plot(hist)

#- convert complex char to binaries

hist$tip.label
hist$states
tb <- convert_to_binary_chars(hist$states)
tb
tb[,1]

# inference

Qsym <- initQ(c(0,1), c(1), diag.as = NA)
Qas <- initQ(c(0,1), c(1,2), diag.as = NA)

# Inference
trait <- cbind(rownames(tb), tb[,5])
res.sym <- rayDISC(tree, trait, rate.mat=Qsym, node.states="none",  model="ARD", root.p="maddfitz")
res.as <- rayDISC(tree, trait, rate.mat=Qas, node.states="none",  model="ARD", root.p="maddfitz")

res.sym$AIC
res.as$AIC

