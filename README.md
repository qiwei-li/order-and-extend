# order-and-extend
```r
#source all functions
source("OE.R")

#the fully known matrix, the example here is a picture
true_matrix = readPNG('image.png')[,,1] 

#create an incomplete matrix, visibility is 30%
visible_matrix = OE_visible(true_matrix = true_matrix, visibility = 0.3) 

#set addtional query budge to 500
query_budge = 500

#set the assumed rank to 2
r=2

#create the order
PI <- OE_order(visible_matrix = visible_matrix, r, adjustment = FALSE)

#run the algorithm
answer <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.1)

#result
answer[[1]] is the error
answer[[2]] is the recovered matrix
answer[[3]] is the positions that the algorithm actively queried