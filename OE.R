# Order and Extend code
# Qiwei Li, UC Davis
# the algorithm is invented in the paper : Matrix Competion With Queries
# this is my own implementation in R

library("Matrix")
library("MASS")
myplot = function(visible_matrix){
  visible_matrix[!is.na(visible_matrix)] = visible_matrix[!is.na(visible_matrix)] + abs(min(visible_matrix, na.rm = TRUE))
  image(t(visible_matrix[nrow(visible_matrix):1, ]), col=grey(seq(0,1,length=256)), axes=FALSE)
}

simulate_matrix = function(r, n1, n2){
  X = matrix(runif(r*n1), nrow=n1)
  Y = matrix(runif(r*n2), nrow=r)
  M = X%*%Y
  X_name = paste0("X", c(1:nrow(M)))
  Y_name = paste0("Y", c(1:ncol(M)))
  rownames(M) <- X_name
  colnames(M) <- Y_name
  return(list(X=X, Y=Y, M=M))
}

OE_stabilize = function(PIi, q1, candidate, visible_matrix, A, t, X, Y, condition_limit){
  #browser()
  if(length(candidate)==0){
    return(NULL)
  }
  candidate = candidate[!candidate %in% q1]
  if(length(candidate)==0){
    return(NULL)
  }
  candidate_condition = rep(9999, length(candidate))
  for(i in 1:length(candidate)){
    if(PIi[3] == "X"){
      pool = c(visible_matrix[as.numeric(PIi[4]), ], visible_matrix[ ,candidate[i]])
      pool = pool[!is.na(pool)]
      if(length(pool)==1){
        tao = pool
      } else {
        tao = sample(x = pool, size = 1)
      }
      A2 = rbind(A, Y[, candidate[i]])
    } else {
      pool = c(visible_matrix[candidate[i], ], visible_matrix[ ,as.numeric(PIi[4])])
      pool = pool[!is.na(pool)]
      if(length(pool)==1){
        tao = pool
      } else {
        tao = sample(x = pool, size = 1)
      }
      A2 = rbind(A, X[candidate[i],])
    }
    t2 = rbind(t, tao)
    s2 = ginv(t(A2) %*% A2) %*% t(A2) %*% t2
    candidate_condition[i] = base::norm(ginv(A2), "2") * base::norm(t2, "2") / base::norm(s2, "2")
  }
  if(min(candidate_condition) < condition_limit){
    numberOfTry <<- 0
    return(candidate[which.min(candidate_condition)])
  } else {
    numberOfTry <<- numberOfTry + 1
    return(NULL)
  }
}

OE_extend = function(PI,  visible_matrix, true_matrix, r, query_budge, condition_limit){
  #browser()
  numberOfTry <<- 0
  orig_visible = visible_matrix
  XY = OE_initialize(PI, r, nrow(true_matrix), ncol(true_matrix))
  X=XY$X
  Y=XY$Y
  PI = PI[-c(1:(2*r)), ]
  while(nrow(PI)>5 & numberOfTry < nrow(PI)){
    print(paste0("Budge: ", query_budge, " Try: ", numberOfTry, " PI: ", nrow(PI)))
    if(PI[1,3]=="X"){
      if(sum(is.na(X[as.numeric(PI[1,4]), ]))==0){
        PI = PI[c(2:nrow(PI), 1), ]
        next
      }
      visible_T = which(!is.na(visible_matrix[as.numeric(PI[1,4]), ]))
      visible_Y = which(colSums(is.na(Y))==0)
      available_equation = intersect(as.numeric(visible_T), as.numeric(visible_Y))
      candidate = visible_Y[which(!visible_Y %in% available_equation)]
      
      # we need to complete the system first, then stable
      if(length(available_equation)<r){
        # not enough equation
        if(length(candidate)+length(available_equation)<r){
          q = "not enough candidate"
        # just enough
        } else if(length(candidate)+length(available_equation)==r){
          if(length(candidate)==1){
            if(query_budge<1){break}
            q = candidate
            query_budge = query_budge-1
          } else {
            if(query_budge<r-length(available_equation)){break}
            q = sample(x = candidate, size= r-length(available_equation))
            query_budge = query_budge-(r-length(available_equation))
          }
          A = t(Y[ , c(available_equation,q)]) # a+q by r
          visible_matrix[as.numeric(PI[1,4]), q] = true_matrix[as.numeric(PI[1,4]), q]
          t = t(t(visible_matrix[as.numeric(PI[1,4]), c(available_equation,q)])) # a+q by 1
          s = t(ginv(A) %*% t) # ginv(A): r by a+q, s = 1 by r
          condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
          if(condition < condition_limit){
            X[as.numeric(PI[1,4]), ] = s
            PI = PI[-1, ]
          } else {
            PI = PI[c(2:nrow(PI), 1), ]
            next #used up the last candidate, still unstable, move on
          }
        # have more candidate
        } else {
          if(query_budge<r-length(available_equation)){break}
          q = sample(x = candidate, size= r-length(available_equation))
          query_budge = query_budge-(r-length(available_equation))
          A = t(Y[ , c(available_equation,q)]) # a+q by r
          visible_matrix[as.numeric(PI[1,4]), q] = true_matrix[as.numeric(PI[1,4]), q]
          t = t(t(visible_matrix[as.numeric(PI[1,4]), c(available_equation,q)])) # a+q by 1
          s = t(ginv(A) %*% t) # ginv(A): r by a+q, s = 1 by r
          condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
          if(condition < condition_limit){
            X[as.numeric(PI[1,4]), ] = s
            PI = PI[-1, ]
          } else {
            a = OE_stabilize(PI[1,], q, candidate, visible_matrix, A, t, X, Y, condition_limit)
            if(is.null(a)){
              PI = PI[c(2:nrow(PI), 1), ]
              next
            }
            if(query_budge<1){break}
            visible_matrix[as.numeric(PI[1,4]), a] = true_matrix[as.numeric(PI[1,4]), a]
            query_budge = query_budge-1
            tao = true_matrix[as.numeric(PI[1,4]), a]
            A2 = rbind(A, Y[,a])
            t2 = rbind(t, tao)
            s2 = ginv(t(A2) %*% A2) %*% t(A2) %*% t2
            X[as.numeric(PI[1,4]), ] = s2
            PI =PI[-1, ]
          }
        }
      # already complete, now stable
      } else {
        A = t(Y[ , c(available_equation)])
        t = t(t(visible_matrix[as.numeric(PI[1,4]), c(available_equation)]))
        s = ginv(t(A) %*% A) %*% t(A) %*% t
        condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
        if(condition < condition_limit){
          X[as.numeric(PI[1,4]), ] = s
          PI = PI[-1, ]
        } else {
          q = NULL
          a = OE_stabilize(PI[1,], q, candidate, visible_matrix, A, t, X, Y, condition_limit)
          if(is.null(a)){
            PI = PI[c(2:nrow(PI), 1), ]
            next
          }
          if(query_budge<1){break}
          visible_matrix[as.numeric(PI[1,4]), a] = true_matrix[as.numeric(PI[1,4]), a]
          query_budge = query_budge-1
          tao = true_matrix[as.numeric(PI[1,4]), a]
          A2 = rbind(A, Y[,a])
          t2 = rbind(t, tao)
          s2 = ginv(t(A2) %*% A2) %*% t(A2) %*% t2
          X[as.numeric(PI[1,4]), ] = s2
          PI = PI[-1, ]
        }
      }
    }
    
    if(PI[1,3]=="Y"){
      if(sum(is.na(Y[ ,as.numeric(PI[1,4])]))==0){
        PI = PI[c(2:nrow(PI), 1), ]
        next
      } 
      visible_T = which(!is.na(visible_matrix[ ,as.numeric(PI[1,4])]))
      visible_X = which(rowSums(is.na(X))==0)
      available_equation = intersect(visible_T, visible_X)
      candidate = visible_X[which(!visible_X %in% available_equation)]

      # we need to complete the system first, then stable
      if(length(available_equation)<r){
        # not enough equation
        if(length(candidate)+length(available_equation)<r){
          q = "no enough candidate"
        }
        # just enough
        else if(length(candidate)+length(available_equation)==r){
          if(length(candidate)==1){
            if(query_budge<1){break}
            q = candidate
            query_budge = query_budge-1
          } else {
            if(query_budge<r-length(available_equation)){break}
            q = sample(x = candidate, size= r-length(available_equation))
            query_budge = query_budge-(r-length(available_equation))
          }
          A = X[c(available_equation,q), ] # a+q by r
          visible_matrix[q, as.numeric(PI[1,4])] = true_matrix[q, as.numeric(PI[1,4])]
          t = t(t(visible_matrix[c(available_equation,q), as.numeric(PI[1,4])])) #a+q by 1
          s = ginv(A) %*% t
          condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
          if(condition < condition_limit){
            Y[, as.numeric(PI[1,4])] = s
            PI = PI[-1, ]
          } else {
            PI = PI[c(2:nrow(PI), 1), ]
            next #used up the last candidate, still unstable, move on
          }
        # have more candidate 
        } else {
          if(query_budge<r-length(available_equation)){break}
          q = sample(x = candidate, size= r-length(available_equation))
          query_budge = query_budge-(r-length(available_equation))
          A = X[c(available_equation,q), ] # a+q by r
          visible_matrix[q, as.numeric(PI[1,4])] = true_matrix[q, as.numeric(PI[1,4])]
          t = t(t(visible_matrix[c(available_equation,q), as.numeric(PI[1,4])])) #a+q by 1
          s = ginv(A) %*% t
          condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
          if(condition < condition_limit){
            Y[, as.numeric(PI[1,4])] = s
            PI = PI[-1, ]
          } else {
            a = OE_stabilize(PI[1,], q, candidate, visible_matrix, A, t, X, Y, condition_limit)
            if(is.null(a)){
              PI = PI[c(2:nrow(PI), 1), ]
              next
            }
            if(query_budge<1){break}
            visible_matrix[a, as.numeric(PI[1,4])] = true_matrix[a, as.numeric(PI[1,4])]
            query_budge=query_budge-1
            tao = true_matrix[a, as.numeric(PI[1,4])]
            A2 = rbind(A, X[a, ])
            t2 = rbind(t, tao)
            s2 = ginv(t(A2) %*% A2) %*% t(A2) %*% t2
            Y[ ,as.numeric(PI[1,4])] = s2
            PI = PI[-1, ]
          }
        }
      # already complete, now stable
      } else {
        A = X[available_equation, ] # a+q by r
        t = t(t(visible_matrix[c(available_equation), as.numeric(PI[1,4])])) 
        s = ginv(t(A) %*% A) %*% t(A) %*% t
        condition = base::norm(ginv(A), "2") * base::norm(t, "2") / base::norm(s, "2")
        if(condition < condition_limit){
          Y[, as.numeric(PI[1,4])] = s
          PI = PI[-1, ]
        } else {
          q=NULL
          a = OE_stabilize(PI[1,], q, candidate, visible_matrix, A, t, X, Y, condition_limit)
          if(is.null(a)){
            PI = PI[c(2:nrow(PI), 1), ]
            next
          }
          if(query_budge<1){break}
          visible_matrix[a, as.numeric(PI[1,4])] = true_matrix[a, as.numeric(PI[1,4])]
          query_budge=query_budge-1
          tao = true_matrix[a, as.numeric(PI[1,4])]
          A2 = rbind(A, X[a, ])
          t2 = rbind(t, tao)
          s2 = ginv(t(A2) %*% A2) %*% t(A2) %*% t2
          Y[ ,as.numeric(PI[1,4])] = s2
          PI = PI[-1, ]
        } #end of stable
      } #end of more eq
    } #end of Y
  } #end of for
  
  recover = visible_matrix
  background = X%*%Y
  recover[is.na(visible_matrix)] = background[is.na(visible_matrix)]
  
  recover2 = recover
  recover2[is.na(recover2)]=0
  
  query = visible_matrix
  query[!is.na(orig_visible)] = NA
  return(list( error = base::norm((true_matrix - recover2), type = "F")/base::norm(true_matrix, type = "F"),
              recover = recover,
              query = query))
}

OE_visible = function(true_matrix, visibility){
  visible_set = sample(x=c(1:length(true_matrix)), size=round(length(true_matrix)*visibility))
  visible_matrix = true_matrix
  visible_matrix[-visible_set] = NA
  return(visible_matrix)
}

OE_initialize = function(PI, r, n1, n2){
  X = matrix(NA, nrow=n1, ncol=r)
  Y = matrix(NA, nrow=r, ncol=n2)
  X[as.numeric(PI[which(PI[,3] =="X")[1:r], 4]), ] = diag(r)
  Y[ ,as.numeric(PI[which(PI[,3] =="Y")[1:r], 4])] = diag(r)
  return(list(X=X, Y=Y))
}

OE_order = function(visible_matrix, r, adjustment=FALSE){
  X_name = paste0("X", c(1:nrow(visible_matrix)))
  Y_name = paste0("Y", c(1:ncol(visible_matrix)))
  rownames(visible_matrix) <- X_name
  colnames(visible_matrix) <- Y_name
  adjacency_matrix = !is.na(visible_matrix)
  indegree = cbind(c(rownames(adjacency_matrix), colnames(adjacency_matrix)), c(rowSums(adjacency_matrix), colSums(adjacency_matrix)))
  indegree = cbind(indegree, gsub("\\d+", '', indegree[,1]), as.numeric(gsub('([A-Z])', '', indegree[,1])))
  first_order = rep(NA, nrow(visible_matrix)+ncol(visible_matrix))
  for(i in length(first_order):1){
    if(0 %in% dim(adjacency_matrix)){
      if(is.null(colnames(adjacency_matrix))){
        first_order[1:length(rownames(adjacency_matrix))] = rownames(adjacency_matrix)
      } else {
        first_order[1:length(colnames(adjacency_matrix))] = colnames(adjacency_matrix)
      }
    }
    X_indegree = rowSums(adjacency_matrix)
    Y_indegree = colSums(adjacency_matrix)
    least = (c(rownames(adjacency_matrix), colnames(adjacency_matrix))[order(c(X_indegree, Y_indegree), decreasing=TRUE)])[i]
    first_order[i] = least
    if(length(grep("X", least))){
      adjacency_matrix = adjacency_matrix[-which(rownames(adjacency_matrix)==least), , drop=FALSE]
    } else {
      adjacency_matrix = adjacency_matrix[, -which(colnames(adjacency_matrix)==least), drop=FALSE]
    }
  }
  first_indegree = indegree[indegree[ ,1][first_order], ]
  if(adjustment==FALSE){
    return(first_indegree)
  } else {
    adjacency_matrix = !is.na(visible_matrix)
    indegree = indegree[indegree[ ,1][first_order], ]
    indegree = cbind(indegree, gsub("\\d+", '', indegree[,1]), as.numeric(gsub('([A-Z])', '', indegree[,1])))
    linked_list = matrix(NA, nrow=2, ncol=length(first_order))
    linked_list[1,]=c("start", first_order[1:(length(first_order)-1)])
    linked_list[2,]=first_order
    for(i in seq_along(first_order)){
      if(indegree[i,2]==0){
        next
      }
      else if(indegree[i,2]>r){
        #it is repositioned to appear immediately after 
        #its neighbor v with the the r-th smallest π(v).
        if(indegree[i,3]=="X"){
          its_neighbor = names(which(adjacency_matrix[as.numeric(indegree[i,4]), ]))
        } else {
          its_neighbor = names(which(adjacency_matrix[ ,as.numeric(indegree[i,4])]))
        }
        if(length(its_neighbor)<r){
          rth_smallest_neighbor = its_neighbor[which.min(sapply(its_neighbor, function(j) which(indegree[,1]==j)))]
        } else {
          rth_smallest_neighbor = its_neighbor[order(sapply(its_neighbor, function(j) which(indegree[,1]==j)), decreasing=TRUE)[r]]
        }
        linked_list = ll_move(linked_list, indegree[i,1], rth_smallest_neighbor)
      } else {
        #it is repositioned to appear immediately after 
        #the neighbor v with the largest π(v).
        if(indegree[i,3]=="X"){
          its_neighbor = names(which(adjacency_matrix[as.numeric(indegree[i,4]), ]))
        } else {
          its_neighbor = names(which(adjacency_matrix[ ,as.numeric(indegree[i,4])]))
        }
        largest_neighbor=its_neighbor[which.min(sapply(its_neighbor, function(j) which(indegree[,1]==j)))]
        linked_list = ll_move(linked_list, indegree[i,1], largest_neighbor)
      }
    }
    second_order = ll_traversal(linked_list)
    second_indegree = indegree[indegree[ ,1][second_order], ]
    return(second_indegree)
  }
}

ll_move = function(linked_list, current, destination){
  # remove current
  current_from = linked_list[1, which(linked_list[2,] == current)]
  current_to = linked_list[2, which(linked_list[1,] == current)]
  linked_list[1, which(linked_list[2,] == current)] = NA
  linked_list[1, which(linked_list[2,] == current_to)] = current_from
  # add "current" after "destination"
  destination_from = linked_list[1, which(linked_list[1,] == destination)]
  destination_to = linked_list[2, which(linked_list[1,] == destination)]
  linked_list[1, which(linked_list[2,] == current)] = destination
  linked_list[1, which(linked_list[2,] == destination_to)] = current
  return(linked_list)
}

ll_traversal = function(linked_list){
  start = which(linked_list[1, ]=="start")
  second_order = rep(NA, ncol(linked_list))
  for(i in 1:ncol(linked_list)){
    second_order[i] = linked_list[2, start]
    start = which(linked_list[1, ]==linked_list[2, start])
  }
  return(second_order)
}