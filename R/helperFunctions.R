# v1.12

require(igraph)
require(arules)


# converts a vector of names into a dataframe
# with original names and corresponding codes A, B, C, D ... , Z, AA, AB, AC, ...
alphaColNames <- function(names) {
  df <- data.frame(names = names, code = NA)
  n <- nrow(df)
  if (n < length(LETTERS)) df$code <- LETTERS[1:n] else 
    if (n > length(LETTERS)*length(LETTERS)) {
      df$code <- paste0("A", as.character(1:n))
    } else {
      df$code[1:length(LETTERS)] <- LETTERS
      n <- n - length(LETTERS)
      ind <- length(LETTERS) + 1
      i <- 1
      while (n > length(LETTERS)) {
        df$code[ind:(ind+length(LETTERS)-1)] <- paste0(LETTERS[i], LETTERS)
        ind <- ind + length(LETTERS)
        i <- i + 1
        n <- n - length(LETTERS)
      }
      df$code[ind:(ind+n-1)] <- paste0(LETTERS[i], LETTERS)[1:n] 
    }
  
  df
  
}

# checks whether a collection is a subset of another collection
# collections are stored as strings, e.g. {A, B, C}

is.subset.of <- function(left, right) {
  left <- substr(left, 2, nchar(left)-1)
  right <- substr(right, 2, nchar(right)-1)
  left <- str_split(left, ",")[[1]]
  right <- str_split(right, ",")[[1]]
  all(left %in% right)
  
}


is.superset.of <- function(left, right) {
  left <- substr(left, 2, nchar(left)-1)
  right <- substr(right, 2, nchar(right)-1)
  left <- str_split(left, ",")[[1]]
  right <- str_split(right, ",")[[1]]
  all(right %in% left)
  
}


intersection.of <- function(left, right) {
  #print("INTERSECTION")
  #print(left)
  #print(right)
  left <- substr(left, 2, nchar(left)-1)
  right <- substr(right, 2, nchar(right)-1)
  left <- str_split(left, ",")[[1]]
  right <- str_split(right, ",")[[1]]
  intersect(left, right)
}

setdiff.of <- function(left, right) {
  #print("SETDIFF")
  #print(left)
  #print(right)
  left <- substr(left, 2, nchar(left)-1)
  right <- substr(right, 2, nchar(right)-1)
  left <- str_split(left, ",")[[1]]
  right <- str_split(right, ",")[[1]]
  res <- setdiff(left, right)
  if (length(res) > 1) res <- paste(res, collapse=",")
  if (length(res) > 0) res <- paste0("{",res,"}")
  res
}



# This are reimplementations of the top down strategy for bulding dendrograms
# out of frequent itemsets collection
# Original implementations are integrated in Python application CollectionWidget


# we are sending in frequent itemsets
# it is expected that quality measures "noItems" and "absSup"
# have already been added
treeTDstrategy <- function(itemsets, root = NA, verbose = F) {
  df <- DATAFRAME(itemsets)
  df$items <- as.character(df$items)
  
  # "bound" parameter tells us whether the collection has been assigned to the tree
  # or if it was permanently discarded
  df$bound <- "UNKNOWN"
  df$bound <- factor(df$bound, levels = c("BOUND", "UNKNOWN", "DISCARDED"), ordered = T)
  df$connectedTo <- NA
  
  # let's rearrange collections so the largest ones are on top
  # same sized are arranged based on absolute support
  # in the case of further tie, we sort them alphabetically
  df <- arrange(df, desc(noItems), desc(absSup), items)
  
  # levels correspond to the size of collections
  # we need to go from largest collections towards the smallest
  currentLevel <- max(df$noItems)
  
  # if we are on the lower levels, re-arrange the collections based on boundness
  # we will first have bound collections, then unknown, then discarded)
  while(currentLevel != 0) {
    
    if (currentLevel != max(df$noItems) - 1) df <- arrange(df, bound)
    
    # let's focus first only on the collections on the current level
    curLevItemsets <- df[df$noItems == currentLevel, "items"]
    
    for (i in 1:length(curLevItemsets)) {
      it <- curLevItemsets[i]
      if (!is.na(root) && !is.subset.of(it, root)) df[df$items == it, "bound"] <- "DISCARDED"
      b <- df[df$items == it, "bound"]
      #if (verbose) print(df)
      
      if (b!="DISCARDED") {
        
        # check out collections on the same level
        # if they share elements with the one we are currently watching
        # discard them
        
        if (i < length(curLevItemsets)) {
          remainder <- curLevItemsets[(i+1):length(curLevItemsets)]
          forDiscard <- remainder [ 
                sapply(sapply(remainder, intersection.of, it), length) > 0]
          df[df$items %in% forDiscard & df$bound == "UNKNOWN", "bound"] <- "DISCARDED"
          
        }
        
      
        # now that we have discarded all collection on the same level with common elements
        # we must move towards the bottom
          
        lowerLevel <- currentLevel - 1
        candidates <- character(0)
        lowLevItemsets <- df[df$noItems == lowerLevel, "items"]
        it.cpy <- it  # we will gradually remove items from the copy
        
        # we now move through lower levels
        # if a collection is still unbound and a subset of it.cpy
        # we will set it as a candidate for binding to the top collection
        # and remove the elements from the copy (since they are now "spent")
        # if it's not a subset, but has common elements, we discard it
        
        
        while (lowerLevel != 0) {
          for (j in 1:length(lowLevItemsets)) {
            lowIt <- lowLevItemsets[j]
            if (verbose) print(paste0("Radim sa podretkom ", lowIt))
            if (df[df$items == lowIt, "bound"] == "UNKNOWN") {
              if (is.subset.of(lowIt, it.cpy)) {
                candidates <- c(candidates, lowIt)
                it.cpy <- setdiff.of(it.cpy, lowIt)
              } else if ((is.subset.of(lowIt, it) == F) && length(intersection.of(lowIt, it)) > 0 ) {
                df[df$items == lowIt, "bound"] <- "DISCARDED"
              }
            }
            if (length(it.cpy)== 0) break # all elements spent
          }
          if (length(it.cpy)== 0) break
          lowerLevel <- lowerLevel - 1
          lowLevItemsets <- df[df$noItems == lowerLevel, "items"]
        }
        
        # once all elements are spent we bind the candidates
        # as well as the top collection
        if (length(it.cpy)==0) {
            df[df$items == it, "bound"] <- "BOUND"
            df[df$items %in% candidates, "bound"] <- "BOUND"
            df[df$items %in% candidates, "connectedTo"] <- it
          } else 
            df[df$items == it & df$bound != "BOUND", "bound"] <- "DISCARDED"
          
          
        }
                              
      }
    currentLevel <- currentLevel - 1
    }
  
  # we now convert solution to the graph
  df.bound <- df[df$bound == "BOUND", ]
  df.bound <- arrange(df.bound, desc(noItems), desc(absSup), items)
  g <- graph.empty(nrow(df.bound), directed = TRUE)
  V(g)$name <- df.bound$items
  #V(g)$name <- substring(V(g)$name, 2, nchar(V(g)$name)-1)
  if (verbose) print(df)
  for (i in 1:nrow(df.bound)) {
    if (!is.na(df.bound[i, "connectedTo"]))
      g <- add.edges(g, edges = c(df.bound[i, "items"], df.bound[i, "connectedTo" ]))
  }
  V(g)$fullName <- df.bound$fullName
  E(g)$weight <- 1
  g
}


# reimplementation of the BU strategy
treeBUstrategy <- function(itemsets, root = NA, verbose = F) {
  df <- DATAFRAME(itemsets)
  df$items <- as.character(df$items)
  
  # "bound" parameter tells us whether the collection has been assigned to the tree
  # or if it was permanently discarded
  df$bound <- "UNKNOWN"
  df$bound <- factor(df$bound, levels = c("BOUND", "UNKNOWN", "DISCARDED"), ordered = T)
  df$connectedTo <- NA
  
  # let's rearrange collections so the largest ones are on top
  # same sized are arranged based on absolute support
  # in the case of further tie, we sort them alphabetically
  df <- arrange(df, desc(noItems), desc(absSup), items)
  
  # we start with two-element collections
  # (singleton colections are all included)
  currentLevel <- 2
  
  while (currentLevel != (max(df$noItems) + 1)) {
    curLevItemsets <- df[df$noItems == currentLevel, "items"]
    if (verbose) {
      print("--------------------------------")
      print("Working on the following collections:")
      print(curLevItemsets)
    }
    
    
    for (i in 1:length(curLevItemsets)) {
      candidates <- character(0)
      
      it <- curLevItemsets[i]
      if (!is.na(root) && !is.subset.of(it, root)) {
        df[df$items == it, "bound"] <- "DISCARDED"
        next
      }
      
      
      it.cpy <- it
      underLevel <- currentLevel - 1
      if (verbose) {
        print(paste0("currentlevel:", as.character(currentLevel)))
        print(paste0("underlevel:", as.character(underLevel)))
        
      }
      
      while (underLevel != 0) {
        # checking out collections on lower levels
        undLevItemsets <- df[df$noItems == underLevel, "items"]
        if (verbose) {
          print("+++++++++++++++++++++++++++++++")
          print("Working on the following undercollections:")
          print(undLevItemsets)
        }
        for (j in 1:length(undLevItemsets)) {
          lowIt <- undLevItemsets[j]
          if ((is.na(df[df$items == lowIt, "connectedTo"])) && is.subset.of(lowIt, it.cpy) &&
              (df[df$items == lowIt, "bound"] != "DISCARDED")) {
            if (verbose) {
              print(paste0("Collection ", lowIt, " is a subset of ", it.cpy))
              print("Adding it to candidates!")
            }
            candidates <- c(candidates, lowIt)
           
            it.cpy <- setdiff.of(it.cpy, lowIt)
            if (verbose) {
              print("Candidates are now:") 
              print(candidates)
              print(paste0("And it.cpy is:", it.cpy))
            }
            
          }
          if (verbose && length(it.cpy)== 0) 
            print("it.cpy is empty! Breaking out!")
          if (length(it.cpy)== 0) break
        }
        if (length(it.cpy)== 0) break
        underLevel <- underLevel - 1
      }
     
      
      # if it.cpy didn't spend all the elements, discard it
      # else bind candidates and it together
      if (length(it.cpy) > 0) {
   
        df[df$items == it, "bound"] <- "DISCARDED"
        if (verbose) {
          print(paste0("Since it.cpy was not empty, I am discarding ", it))
          print("Our data frame now looks like this:")
          print(df)
        }
        
      } else {
        df[df$items == it, "bound"] <- "BOUND"
        df[df$items %in% candidates, "bound"] <- "BOUND"
        df[df$items %in% candidates, "connectedTo"] <- it
        
        if (verbose) {
          print("it.cpy is empty! We are bounding the candidates.")
          print("Our data frame now looks like this:")
          print(df)
        }
      }
    
      
    }
    
    currentLevel <- currentLevel + 1
    
  }
  
  df[df$bound == "UNKNOWN", "bound"] <- "DISCARDED"

  # we now convert solution to the graph
  df.bound <- df[df$bound == "BOUND", ]
  df.bound <- arrange(df.bound, desc(noItems), desc(absSup), items)
  if (verbose) print(df)
  g <- graph.empty(nrow(df.bound), directed = TRUE)
  V(g)$name <- df.bound$items
  #V(g)$name <- substring(V(g)$name, 2, nchar(V(g)$name)-1)
  for (i in 1:nrow(df.bound)) {
    if (!is.na(df.bound[i, "connectedTo"]))
      g <- add.edges(g, edges = c(df.bound[i, "items"], df.bound[i, "connectedTo" ]))
  }
  V(g)$fullName <- df.bound$fullName
  E(g)$weight <- 1
  g
}


# this is a variation of a function from 'arules' package
# the procedure is slightly tweaked (so it doesn't break down
# when multiple smaller trees are built)
# and it returns an igraph now
msTreeKruskalCustom <- function (itemsets, root = NA, minSize = 1, removeLoneSingletons = TRUE) {
  
  # if we want to reduce the size
  t <- itemsets[quality(itemsets)$noItems >= minSize]
  
  #if we have a root
  if (!is.na(root)) {
    t <- itemsets[sapply(as.character(DATAFRAME(itemsets)$items),  is.subset.of, root), ]
    if (length(t)==0) t <- itemsets  # wrong root, silently go forward
  }
  
  if (minSize == 1 && removeLoneSingletons) {
    
    # let's remove singletons who don't have supersets
    # WARNING: awkward code incoming
    tt <- DATAFRAME(itemsets)
    tt$items <- as.character(tt$items)
    tt1 <- tt[tt$noItems == 1, ]
    
    # basically, what we have here should be a nested for loop
    # for each singleton collection I need to iterate over all the non-singletons
    # and check if one is a superset
    # since this is R, I will use two sapplys
    # the first one takes a list of singletons
    # and applies a function which checks for supersettage
    # this function contains a second sapply which uses a non-singleton list
    # and the original singleton as a parameter
    # I call this "sapplyception"
    # everything is ultimately wrapped in a data frame with a singleton name
    # and whether it should be kept
    tt1 <- data.frame(items = tt1$items, keepMeIHaveASuperset = sapply(tt1[, 1], 
                          function(x) any(sapply(tt[tt$noItems > 1,1], is.superset.of, x))),
                      stringsAsFactors = F)
    tt <- left_join(tt, tt1, by = "items")$keepMeIHaveASuperset
    tt[is.na(tt)] <- T
    
    t <- itemsets[tt]
    
  }
  
  # we create an upper triangle matrix where we fill only
  # the cells where one itemset is a subset of another
  #
  # the value we put in is the absolute difference between 
  # absolute support values
  m <- rep(NA, length(t)^2)
  m <- matrix(m, ncol = length(t))  # square matrix of NAs
  res <- DATAFRAME(t)
  rownames(m) <- res[,1]            # collections are rows and columns
  colnames(m) <- res[,1]
  
  # loops in R are bad! (but sometimes we like to have them anyway)
  
  for (i in 1:(length(t)-1))
    for (j in (i+1):length(t))
      if (as.logical(is.subset(t[i], t[j])) || as.logical(is.subset(t[j], t[i]))) {
        m[i,j] <- abs(quality(t[i])$absSup - quality(t[j])$absSup)
        #print(paste0("i je ", i, " a j je ", j, " a m[i,j] je ", m[i,j]))
      }
  
  # now m is still a matrix of NAs, except on cells where the column collection
  # is a subset of row collection
  # in that case the cell contains the difference in absolute support between
  # the subset and the superset
  
  #let's build a graph, nodes are collections
  nodes <- res[,1]
  arcs <- NULL
  # again, this should preferably be in a more R-friendly format
  # (and some day, it will be :) )
  for (i in 1:(length(t)-1))
    for (j in (i+1):length(t))
      if (!is.na(m[i,j])) {
        # for each cell which isn't NA we add an arc 
        if (is.null(arcs)) {
          arcs <- c(i, j, m[i,j])
        } else arcs <- rbind(arcs, c(i, j, m[i,j]))
      }
  
  # we order the arc based on the support difference
  # most of the code right below is from the original Kruskal function
  arcs <- matrix(arcs[order(arcs[, 3]), ], ncol = 3)
  components <- matrix(c(nodes, nodes), ncol = 2)
  tree.arcs <- matrix(ncol = 3)[-1, ]
  stages <- 0
  stages.arcs <- c()
  i <- 1
  while (nrow(tree.arcs) < length(nodes) - 1 && i <= nrow(arcs)) {
    min.arc <- arcs[i, ]
    iComp <- components[components[, 1] == min.arc[1], 2]
    jComp <- components[components[, 1] == min.arc[2], 2]
    if (iComp != jComp) {
      tree.arcs <- rbind(tree.arcs, min.arc)
      components[components[, 2] == jComp, 2] <- iComp
    }
    stages <- stages + 1
    stages.arcs <- c(stages.arcs, rep(stages, nrow(tree.arcs) - 
                                        length(stages.arcs)))
    i <- i + 1
  }
  
  # let's delet the row names and call the columns ept1, ept2 and weight
  # simply because the original algorithm did it like this (I think)
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
  rownames(tree.arcs) <- NULL
  
  # original return of the function
  # output<- list(tree.nodes = nodes, tree.arcs = tree.arcs, stages = stages, stages.arcs = stages.arcs)
  
  # but we want an Igraph!
  tree.df <- as.data.frame(tree.arcs)
  
  # remember res? It's our original frequent itemsets data frame
  # and the row number of tree.df corresponds with the original row number
  # we can use this to convert the numbers to actual names
  
  # now this code is a little ugly and relies to much on inner joins
  # one day it will become more beautiful
  res$id <- 1:nrow(res)
  tree.df <- inner_join(x = res[, c("items", "id")], y = tree.df, by = c("id" = "ept1"))
  names(tree.df)[1:2] <- c("ept1_name", "ept1_id")
  tree.df <- inner_join(x = res[, c("items", "id")], y = tree.df, by = c("id" = "ept2"))
  names(tree.df)[1:2] <- c("ept2_name", "ept2_id")
  tree.df <- tree.df[, c("ept1_id", "ept1_name", "ept2_id", "ept2_name", "weight")]
  
  
  # let's put larger collections to the left
  for (i in 1:nrow(tree.df)) {
    left <- str_count(tree.df[i, "ept1_name"], ",")
    right <- str_count(tree.df[i, "ept2_name"], ",")
    if (left < right) 
      tree.df[i,] <- tree.df[i, c("ept2_id", "ept2_name", "ept1_id", "ept1_name")]
  }
  
  # time do build a graph
  g <- graph.empty(length(nodes), directed = TRUE)
  V(g)$name <- as.character(nodes)
  V(g)$name <- substring(V(g)$name, 2, nchar(V(g)$name)-1)  # {C,D} -> C,D
  V(g)$fullName <- character(length(V(g)))   # all full names are empty strings
  V(g)$noItems <- str_count(V(g)$name, ",") + 1   # item count can be inferred by the name of the collection
  V(g)$absSup <- res$absSup     # why not add the absolute and relative support?
  V(g)$support <- res$support
  
  # singleton collections get the privilege of having full names
  for (i in 1:length(V(g)))
    if (V(g)[i]$noItems == 1) V(g)[i]$fullName <- 
    as.character(dfColEncodings[dfColEncodings$code == V(g)[i]$name, "names"])
  
  
  # first two columns of tree.arcs are edges
  # since R stores a matrix by columns we transpose the two columns
  # and turn them into a vector, so the edges are defined by a row
  # of numbers, each pair depicting an edge
  g <- add.edges(g, edges = as.vector(t(tree.df[,c("ept1_id", "ept2_id")])))
  E(g)$weight <- tree.df$weight
  g
}