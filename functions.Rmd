#find distance between sequences
distance <- function(seq1, seq2) {
  
  dist <- 0
  for(i in 1:length(seq1)) {

      if(seq1[i] != seq2[i]) {
        dist <- dist + 1
      }
  }
  return(dist)
}

#find distance between all sequences
distanceOfAllPairs<- function(alignmentPatterns, alignmentSubjects){
  
  n <- length(sequences)
  a <- length(alignmentPatterns)
  distVector <- numeric(a)
  
  for(i in 1:a) {
    
    patterns <- strsplit(alignmentPatterns[i], "")[[1]]
    subjects <- strsplit(alignmentSubjects[i], "")[[1]]
      
    distVector[i] <- distance(patterns, subjects)
  }
  
  return(distVector)

}

#create distance matrix
distMatrix <- function(alignmentSubjects, alignmentPatterns, substitutionMatrix, gapOpening, gapExtension){

  distVector <- distanceOfAllPairs(alignmentPatterns,alignmentSubjects)
  n <- length(sequences)
  
  
  distMatrix = matrix(0, nrow = n, ncol = n) #initalize distMatrix with zeros
  dimnames(distMatrix) <- list(sequences, sequences) #label columns and rows according to the corresponding sequence
  distMatrix[lower.tri(distMatrix)] <- distVector #populate lower triangle with distances between sequences
  distMatrix <- distMatrix + t(distMatrix) #sum with own transposition to mirror values across diagonal
  return(distMatrix)
}

#align all sequences with all others
uniqueAlignedCombo <- function(sequences, substitutionMatrix, gapOpening, gapExtension){
  n <- length(sequences)
  alignmentPatterns <- character()
  alignmentSubjects <- character()
  k <- 1
  
  for(i in 1:(n-1)) { 
      for(j in (i + 1):n) { 
          alignmentPatterns[k] <- pairAlign(sequences[i], sequences[j], substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension)[[1]]
          alignmentSubjects[k] <- pairAlign(sequences[i], sequences[j], substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension)[[2]]
          k <- k + 1 
      }
  }
  return(list(alignmentPatterns = alignmentPatterns, 
              alignmentSubjects = alignmentSubjects))
}

#find center star sequence
findCenter <- function(sequences, substitutionMatrix, gapOpening, gapExtension){
  
  #finds pattern-subject pairs in sequences
  pairs <- uniqueAlignedCombo(sequences, substitutionMatrix, gapOpening, gapExtension)
  
  #finds distance matrix for pairs
  dist <- distMatrix(pairs$alignmentSubjects, pairs$alignmentPatterns, substitutionMatrix, gapOpening, gapExtension)
  
  #totals distances in distVec
  distVec <- numeric(ncol(dist))
  for(i in 1:ncol(dist)) {
    distVec[i] <- distVec[i] + sum(dist[,i])
  }
  
  #finds index of minimum value
  c <- which(distVec == min(distVec))
  return(c)

}

#formats inputs for buildMSA
formatForbuildMSA <- function(patternsMat,subjMat, center){ #initializes function with parameters patternsMat, subjMat, and center
  patternsCenter<- patternsMat[ ,center][-center]
  subjCenter<- subjMat[ ,center][-center] 
  
  patternMSA <- strsplit(patternsCenter, "") 
  subjectMSA <- strsplit(subjCenter, "") 
 return(list(patterns = patternMSA, 
             subjects = subjectMSA)) 
}

#creates alignment matrices
createAlignmentMatrices <- function(seq, substitutionMatrix, gapOpening, gapExtension){
  n <- length(seq)
  alignmentPatternsmat <- matrix("", n, n)
  alignmentSubjectsmat <- alignmentPatternsmat
  
  #Generate alignments
  for(i in 1:(n-1)) {
      for(j in (i+1):n) {
          alignment <- pairAlign(seq[i], seq[j], substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension)
          
          alignmentPatternsmat[j, i] <- alignment[[1]]
          alignmentSubjectsmat[i, j] <- alignment[[1]]
          alignmentSubjectsmat[j, i] <- alignment[[2]]
          alignmentPatternsmat[i, j] <- alignment[[2]]
      }
  }
   
   dimnames(alignmentPatternsmat) <- list(seq, seq)
   dimnames(alignmentSubjectsmat) <- list(seq, seq)
  return(list(alignmentPatterns = alignmentPatternsmat, 
              alignmentSubjects = alignmentSubjectsmat))
}

#builds multiple sequence alignment formatting
buildMSA <- function(patterns, subjects, center) {
    MSA <- rbind(patterns[[1]], subjects[[1]])
    for(i in 2:length(patterns)) {
        j = 1 #index in new row
        k = 1 #index in alignment of center sequence to sequence i
        m = 1 #column index of MSA
        maxLength = ncol(MSA) + length(patterns[[i]])
        newRow = character(maxLength)
        while(k <= length(patterns[[i]]) && m <= ncol(MSA)) {
            if(patterns[[i]][k] == MSA[1, m]) {
                newRow[j] <- subjects[[i]][k]
                j <- j + 1
                k <- k + 1
                m <- m + 1
            } else if(MSA[1, m] == "-") {
                newRow[j] <- "-"
                j <- j + 1
                m <- m + 1
            } else if(patterns[[i]][k] == "-") {
                if(m == 1) {
                    MSA <- cbind(rep("-", nrow(MSA)), MSA)
                } else {
                    MSA <- cbind(MSA[, 1:(m-1)], rep("-", nrow(MSA)), MSA[, m:ncol(MSA)])
                }
                newRow[j] <- subjects[[i]][k]
                m <- m + 1
                j <- j + 1
                k <- k + 1
            }
        }
        while(k <= length(patterns[[i]])) {
            MSA <- cbind(MSA, rep("-", nrow(MSA)))
            newRow[j] <- subjects[[i]][k]
            k <- k + 1
            j <- j + 1
        }
        while(m <= ncol(MSA)) {
            newRow[j] <- "-"
            m <- m + 1
            j <- j + 1
        }
        newRow <- newRow[1:(m - 1)]
        MSA <- rbind(MSA, newRow)
    }
    rownames(MSA) <- c("    Center:", paste0("Sequence ", 1:(nrow(MSA)), ":")[-center])
    colnames(MSA) <- 1:ncol(MSA)
    return(MSA)
}

#final function - outputs multiple sequence alignment
centerStar <- function(sequences, substitutionMatrix, gapOpening, gapExtension) {
  
  #finds center
  center <- findCenter(sequences, substitutionMatrix, gapOpening, gapExtension)
  
  #creates pattern + subject matrices
  patMat <- createAlignmentMatrices(sequences, substitutionMatrix, gapOpening, gapExtension)$alignmentPatterns
  subjMat <- createAlignmentMatrices(sequences, substitutionMatrix, gapOpening, gapExtension)$alignmentSubjects
    
  #format
  format <- formatForbuildMSA(patMat, subjMat, center)
  
  #msa
  msa <- buildMSA(format$patterns, format$subjects, center)
  return(msa)
  
  
}
