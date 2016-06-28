# order-and-extend
R implementation of the active matrix completion algorithm in the paper ”Matrix Completion with Queries”

## Example
```r
source("OE.R")

#the fully known matrix and create an incomplete matrix, visibility is 30%
true_matrix = readPNG('image.png')[,,1] 
visible_matrix = OE_visible(true_matrix = true_matrix, visibility = 0.3) 

#set addtional query budge, assumed rank, and run the algorithm
query_budge = 500
r=2
PI <- OE_order(visible_matrix = visible_matrix, r, adjustment = FALSE)
answer <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, 
					condition_limit=1.1)

#result
answer[[1]] is the error
answer[[2]] is the recovered matrix
answer[[3]] is the positions that the algorithm actively queried
```

## Results
![ex](pic1.png)
![ex](pic2.png)

##License
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
