sample_states <- function(total, N){
  sample(c(2:total), size=N)
}

adj2Q <- function(Qr){
  Scaled <- Qr/ rowSums(Qr)
  diag(Scaled) <- -rowSums(Scaled)
  return(Scaled)
}


reorder_strings <- function(strings, position) {
  order(substr(strings, position, position))
}

# Example usage:
strings <- c("000", "001", "010", "011", "100", "101")
position <- 2
reorder_strings(strings, position)

reoder_Q <- function(Q, position) {
  strings <- rownames(Q)
  v <- reorder_strings(strings, position)
  Q[v,v]
}

get_partition_index <- function(Q.reod, position) {
  strings <- rownames(Q.reod)
  x <- substr(strings, position, position)
  which(x==1)[1]
}  


get_lumpQrates <- function(Q.reod, position){
  x1.1=get_partition_index(Q.reod, position)
  x1.2=ncol(Q.reod)
  y1.1=1
  y1.2=x1.1-1
  submatrixR <- Q.reod[c(y1.1:y1.2), c(x1.1:x1.2)]
  x2.1=1
  x2.2=x1.1-1
  y2.1=x1.1
  y2.2=x1.2
  submatrixL <- Q.reod[c(y2.1:y2.2), c(x2.1:x2.2)]
  rR <- sum(submatrixR)/nrow(submatrixR)
  rL <- sum(submatrixL)/nrow(submatrixL)
  return(c(rL, rR))
}


get_lumpQrates_all <- function(Q.reod, nchar=3){
  r <- matrix(NA, nrow=nchar, ncol=2)
  r[1,] <- get_lumpQrates(Q.reod,1)
  for (i in 2:nchar) {
    Qi <- reoder_Q(Q.reod, i)
    r[i,] <- get_lumpQrates(Qi, i)
  }
  r
} 


convert_to_binary_chars <- function(strings) {
  # Split each string into individual characters
  binary_chars <- strsplit(strings, "")
  
  # Determine the number of columns (equal to the number of characters in the strings)
  num_cols <- nchar(strings[[1]])
  
  # Create a matrix to store the binary characters
  binary_table <- matrix(nrow = length(strings), ncol = num_cols)
  
  # Fill in the matrix with the binary characters
  for (i in 1:length(strings)) {
    binary_table[i, ] <- unlist(binary_chars[[i]])
  }
  #binary_table <- cbind(names(strings), binary_table)
  rownames(binary_table) <- names(strings)
  # Return the binary table
  return(binary_table)
}

# Example usage:
strings <- c("00", "01", "10", "11")
binary_table <- convert_to_binary_chars(strings)
print(binary_table)




get_neighbours <- function(Qj, focal.state='111', max.dist){
  states <- rownames(Qj)
  distances <- stringdistmatrix(states, states, method = "hamming")
  focal.state <- which(states==focal.state)
  index=which(distances[focal.state,]==max.dist)
  out <- index
  names(out) <- states[index]
  out
}

get_neighbours_maxdist <- function(Qj, focal.state='111', max.dist){
  nei <- c()
  for (i in 1:max.dist){
     x <- get_neighbours(Qj, focal.state, max.dist=i)
     nei <- c(nei, x)
  }
  # include focal state itself
  states <- rownames(Qj)
  f.inx=which(states==focal.state)
  names(f.inx)=focal.state
  nei <- c(nei, f.inx)
  #
  nei
}



# Function to randomly mutate a binary string for N generations
mutate_binary_string <- function(input_string, generations) {
  mutants <- vector("list", length = generations)
  mutants[[1]] <- input_string
  
  for (gen in 2:generations) {
    chars <- unlist(strsplit(mutants[[gen - 1]], ""))
    mutate_pos <- sample(length(chars), 1)
    chars[mutate_pos] <- as.character(1 - as.numeric(chars[mutate_pos]))
    mutants[[gen]] <- paste(chars, collapse = "")
  }
  
  return(unlist(mutants))
}

# Example usage:
initial_string <- "111"
generations <- 10
mutants <- mutate_binary_string(initial_string, generations)
print(mutants)


# Function to randomly mutate a binary string for N generations, ensuring uniqueness
mutate_unique_binary_string <- function(input_string, generations) {
  mutants <- vector("list", length = generations)
  mutants[[1]] <- input_string
  seen_strings <- c(input_string)
  
  for (gen in 2:generations) {
    repeat {
      chars <- unlist(strsplit(mutants[[gen - 1]], ""))
      mutate_pos <- sample(length(chars), 1)
      chars[mutate_pos] <- as.character(1 - as.numeric(chars[mutate_pos]))
      new_string <- paste(chars, collapse = "")
      
      if (!(new_string %in% seen_strings)) {
        seen_strings <- c(seen_strings, new_string)
        mutants[[gen]] <- new_string
        break
      }
    }
  }
  
  return(unlist(mutants))
}

# Example usage:
# initial_string <- '0000000000'
# generations <- 100
# mutants <- mutate_unique_binary_string(initial_string, generations)
# unique(mutants)
# print(mutants)

