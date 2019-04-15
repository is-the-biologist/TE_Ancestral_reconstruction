library(ape)
library(hash)
library(abind)
library(ggplot2)
library(sandwich)
library(gridExtra)
library(ggrepel)


calc_all_ACE<-function(mei_table){
  
  #This function takes in as input a tabular data structure with presence or absence of an insertion for each of TE call for each strain. A python script will convert VCFs to this format.
  
  
  
  #I want to generate a position x strain x ancestor_node matrix wherein each value of the matrix denotes a change in state from that node regardless if it makes sense phylogenetically -- i can deal with those consequences later this way I should be able to some fairly simple operations summing through the matrix or some such 
  
  eb_tree<- "((((((BL6.NTac:0.04125,BL6.NJ:0.04125):0.01195,(BL6.NCrl:0.04305,BL6N-TyrCBrdCrlCrl:0.04305):0.01015):0.00398,BL6.NHsd:0.05718):0.01931,BL6.ByJ:0.07649):0.00292,(BL6.JBomTac:0.06160,(BL6.J:0.05538,BL6.JEiJ:0.05538):0.00622):0.01781):0.01657,(((BL10.ScSnJ:0.05963,BL10.SnJ:0.05963):0.01815,(BL10.ScCr:0.05382,BL10.ScNHsd:0.05382):0.02397):0.00621, BL10.J:0.08399):0.01199);"
  mouse_tree<- read.tree(text=eb_tree)
  mei_matrix <- as.matrix(mei_table)
  
  
  position_names<- vector()
  
  
  arr_index = 0
  for (insertion in 2:ncol(mei_matrix)){
    ins_trait<-mei_matrix[,insertion]
    
    #In order for ACE to work it needs to consider only traits that are polymorphic between the tips, meaning all 0's or all 1's must be filtered -- all 0's were filtered pre-process but all 1's must be filtered now.
    
    if ( sum(strtoi(ins_trait)) != 14){#Check if vectors sums to 1 if no then proceed
      names(ins_trait)<-mouse_tree$tip.label
      ancestral <- ace(x=ins_trait, phy=mouse_tree, type='discrete', model='SYM')
      #We should change the transition matrix from being symmetrical to one that is specified
      
      
      #The most ancestral node is index 1 
      
      node_ancestor_Probs <- ancestral$lik.anc
      
      
      likelihood <- ancestral$loglik
      
      
      #Now we have a vector of sites and we will iterate through each node and filter by likelihood first
      if (likelihood > -9){
        #If passes our likelihood filter we will check each row of the nodes and see
        ancestor<-ifelse(node_ancestor_Probs[,1] >= .5, 0, 1) #make the value of each node 0 or 1 depending on probability of state
        #get position names
        position_names<-c(position_names, colnames(mei_table)[insertion])
        
        #Now we will check every node against every substrain and note a change in state into a matrix
        #instant position matrix
        pos_matrix <-  matrix(nrow = 27, ncol = 13)
        #Loop for each strain 
        
        for (strain in 1:14){
          strain_trait <- ins_trait[strain]
          #Loop for each ancestor
          #instant array to hold array
          all_changes <- vector()
          for (ancs in 1:13){
            if (ancestor[ancs] == strain_trait){
              #case where ancestor is same state as strain mark as 0
              change_state <- 0
            }
            else if (ancestor[ancs] > strain_trait){
              #case where we have a loss of TE at position
              change_state <- -1
            }
            else if (ancestor[ancs] < strain_trait){
              #case where we have a gain in TE at position
              change_state <- 1
            }
            all_changes <- c(all_changes, change_state)
          }  
          #add state changes to an array 
          pos_matrix[strain,] <- all_changes          
        }
          #We have computed the changes between a given strain and a node, but now we want to compute the differences between a node and any given node
          
          #Iterate through node
        
        for (node_1 in 1:13){
          node_change <- vector()
          for (node_2 in 1:13){
            
            if (ancestor[node_2] == ancestor[node_1]){
              #case where ancestor is same state as strain mark as 0
              change_state <- 0
              
            }
            else if (ancestor[node_2] > ancestor[node_1]){
              #case where we have a loss of TE at position
              change_state <- -1
            }
            else if (ancestor[node_2] < ancestor[node_1]){
              #case where we have a gain in TE at position
              change_state <- 1
            }
            node_change <- c(node_change, change_state)
            
            } 
          
          pos_matrix[node_1+14, ] <- node_change 
          }
          
            
        
          
          
        
        #Array will carry a matrix for each position that contains the changes from each strain to that node
        if (arr_index == 0){
          ancestral_array <- array(data= pos_matrix, dim = c(27, 13, 1))
        }
        else { 
          ancestral_array <- abind(ancestral_array, pos_matrix, along=3)
        }
        arr_index = 1
      }
      
      
    }
    
    
  }
  
  #rownames(ancestor_calls) <- position_names
  #colnames(ancestor_calls) <- branching.times(mouse_tree)
  dimnames(ancestral_array)[[3]] <- position_names
  return(ancestral_array)
  
}


#Function will compute the number of changes from the ancestral nodes of our specification and the 
get_changes<-function(changes_array, node, strains=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)){
  gains_loss_matrix = matrix(data = 0, nrow = 2, ncol = 27)#gains, 1 ; loss, 2
  
  for (position in 1:dim(changes_array)[3]){
    #iterate through each position and calculate the number of gains/losses for each strain 
    
    for (mouse in strains){
      
      state_change <- changes_array[mouse, node, position]
      #check if the state change is a gain, loss or no change:
      if (state_change == 1){
        #gain of state relative to node
        gains_loss_matrix[1, mouse] <- gains_loss_matrix[1, mouse] + 1
        
      }
      else if  (state_change == -1){
        #deletion relative to node
        gains_loss_matrix[2, mouse] <- gains_loss_matrix[2, mouse] + 1
        
      }        
    }
    
  }
  
  return(gains_loss_matrix)
  
  
}


calc_branches <- function(phylogeny, ancestry, branch_lengths, fileNAME){
  ### This function will wrap the get_changes function to get every branch's gains and losses relative to each node and leaf
  #It requires an argument that will match every node/leaf to another node to be calculated
  
  #Tree
  eb_tree<- "((((((BL6.NTac:0.04125,BL6.NJ:0.04125):0.01195,(BL6.NCrl:0.04305,BL6N-TyrCBrdCrlCrl:0.04305):0.01015):0.00398,BL6.NHsd:0.05718):0.01931,BL6.ByJ:0.07649):0.00292,(BL6.JBomTac:0.06160,(BL6.J:0.05538,BL6.JEiJ:0.05538):0.00622):0.01781):0.01657,(((BL10.ScSnJ:0.05963,BL10.SnJ:0.05963):0.01815,(BL10.ScCr:0.05382,BL10.ScNHsd:0.05382):0.02397):0.00621, BL10.J:0.08399):0.01199);"
  mouse_tree<- read.tree(text=eb_tree)
  
  ##################
  BL <- vector()
  gained_branches <- vector()
  lost_branches <- vector()
  for (branches in 1:length(phylogeny)){
    
    #Get all the gains and losses for the two branches connected to that node
    branch_matrix <- get_changes(changes_array = ancestry, node = phylogeny[[branches]][1], strains=c(phylogeny[[branches]][2], phylogeny[[branches]][3]))
    
    gained_branches <- c(gained_branches, branch_matrix[1, phylogeny[[branches]][2]],  branch_matrix[1, phylogeny[[branches]][3]]) #add to gains vector
    lost_branches <- c( lost_branches, branch_matrix[2, phylogeny[[branches]][2]], branch_matrix[2, phylogeny[[branches]][3]])
    BL<-c(BL, mouse_tree$edge.length[branch_lengths[[branches]][1]], mouse_tree$edge.length[branch_lengths[[branches]][2]])
    #now lets add the respective branch lengths to a vector that are respective for those branches that had gains and losses calculated
    
    
    }
 
  
  #Now we have extracted the relevant data we will perform the poisson regression
  TE_df <- data.frame(gains = gained_branches, dels = lost_branches, lengths = BL)
  poiss_gains<-glm(gains ~ lengths, data= TE_df, family='poisson')
  poiss_loss<- glm(dels ~ lengths, data=TE_df, family='poisson')
  
  #calculate p_values of fit
  p_g<-pchisq(poiss_gains$deviance, df=poiss_gains$df.residual, lower.tail = FALSE)
  p_l<-pchisq(poiss_loss$deviance, df=poiss_loss$df.residual, lower.tail = FALSE)
  p_vals <- c(p_g, p_l)
  #calculate predicted gains and deletions based on model
  pred_gains <- predict(poiss_gains, type='response')
  pred_dels <- predict(poiss_loss, type='response')
  TE_df$pred_gains=pred_gains
  TE_df$pred_dels=pred_dels
  
  #Now lets plot and display the regression
  TE_poiss <-ggplot(TE_df) + geom_point( aes(x=lengths, y=gains, color='Insertions')) + xlab(label = 'Branch lengths') + ylab(label = 'TE insertions/deletions') +
    ggtitle(label = paste0(fileNAME, ': insertions/deletions vs. branch length')) + geom_line(aes(x=lengths, y=pred_gains, color='Poisson regression (Insertions)'))  + 
    geom_point(aes(x = lengths, y=dels, color='Deletions')) + geom_line(aes(x=lengths, y=pred_dels, color='Poisson regression (Deletions)')) + 
    theme(text = element_text(size=15), legend.text = element_text(size = 7.5))

  ggsave(paste0('~/Documents/Clark_lab/Poiss_Regressions/', fileNAME, '.png'), plot=TE_poiss)
  #TE_poiss
  
  
  #output p values
  return(p_vals)
}



wrap_regression<- function(){
  #the final wrapper function to calculate the poisson regression on all of our TEs of interest
  
  
  ###################
  #Input in a list of node, connecting node1, connecting node2 
  #remember that when adding the non-leaf nodes as connecting nodes to this list to add 14 to their index as that is the position they are in the matrix (15-27)
  #However the first node object in the vector does not need +14 as it is being taken from the columns of the matrix indexed from 1-13
  phylog<-list(
    c(10, 14, 25),
    c(11, 27, 26),
    c(13, 13, 12),
    c(12, 11, 10),
    c(2, 22, 17),
    c(8, 23, 7),
    c(9, 9, 8),
    c(3, 6, 18),
    c(4, 5, 19),
    c(5, 21, 20),
    c(7, 4, 3),
    c(6, 2, 1)
  )
  #takes the indices of the branch length vector for each of the nodes
  bl<-list(c(26, 19), c(23, 20), c(25, 24), c(22, 21), c(13, 2), c(15, 14), c(17, 16), c(12, 3), c(11, 4), c(8, 5), c(10, 9), c(7, 6)
  ) 
  
  ####################
  p_val_matrix <- matrix(nrow = 5, ncol = 2)
  setwd(dir = '~/Documents/Clark_lab/VCF_outputs/mei_tsv/')
  TE_oi <- c('B1mei.tsv', 'B2mei.tsv', 'L1MdAmei.tsv', 'L1MdGfmei.tsv', 'L1MdTmei.tsv')
  
  for (file in 1:5){
    new_Name<-strsplit(TE_oi[file], 'mei.tsv')
    
    te_table<- read.table(TE_oi[file], sep = '\t', header = TRUE)
    ancestor_table<-calc_all_ACE(mei_table = te_table)
    p<-calc_branches(phylogeny = phylog, ancestry = ancestor_table, branch_lengths = bl, fileNAME = new_Name)
    
    p_val_matrix[file,] <- p
  }
  colnames(p_val_matrix) <- c('gains_deviance', 'deletion_deviance')
  rownames(p_val_matrix) <- TE_oi
  
  write.table(p_val_matrix, file= '~/Documents/Clark_lab/Poiss_Regressions/poiss_pvals.tsv', sep='\t', quote = FALSE)
  
}

wrap_regression()


