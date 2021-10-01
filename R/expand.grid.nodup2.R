#' Combinations of victims without duplication
#'
#' This is based on the function [expand.grid()], but also removes combinations
#' with repeated elements. This function is specifically designed to deal with
#' DVI, and uses DVI data to further cut down on combinations. The function
#' takes sex and age as parameters, and removes impossible combinations based
#' on this information.
#'
#' @param lst A list of vectors.
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param sex A list of the sexes of the victims. (1=male,2=female,0=unknown)
#' @param age A list of age constraints for the victims, these are compared
#' against the structure of the pedigree.
#'
#' @return A data frame. Each row is a possible solution to the disaster victim
#' identification problem.
#'
#' @seealso [expand.grid()]
#'
#' @examples
#' library(stringr) # removed importfrom
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#'
#' am[1] = setAlleles(am[1],"R1", alleles = 0)
#' jointRes = jointDVI(pm,am,missing,verbose = TRUE)
#' 
#' pairss = list(V1=c("*","M1","M2"),V2=c("*","M1","M2"),V3=c("M3","*"))
#' expand.grid.nodup2(pairss,pm,am,age="V2>V1")
#'

#'
#' @export



# NEW_FUNC: function taking age and sex into count.

expand.grid.nodup2 = function(lst,pm=NULL,am=NULL,sex=NULL,age=NULL) {
  
  # copy-paste of old function
  if(is.data.frame(lst))
    stop("Unexpected input: Argument `lst` is a data frame")
  if(!is.list(lst))
    stop("Argument `lst` should be a list")
  
  len = length(lst)
  if(len == 0)
    return(data.frame())
  
  nms = names(lst)
  if(is.null(nms))
    nms = paste0("Var", 1:len)
  
  if(any(lengths(lst) == 0))
    return(data.frame(matrix(nrow = 0, ncol = len), dimnames = list(NULL, nms)))
  
  
  # accounting for sex:
  
  if(!is.null(sex)){
    for(i in 1:len){
      snm = sex[nms[i]]
      if(!is.na(snm)){
        if(snm!=0){
          lsti = lst[[i]]
          rm = c()
          for(j in 1:length(lsti)){
            sli = sex[lsti[j]]
            if(!is.na(sli)){
              if(sli!=0){
                if(sli!=snm){
                  rm = append(rm,i)
                }
              }
            }
          }
          if(!is.null(rm)){
            lst[[i]] = lst[[i]][-rm]
          }
        }
      }
    }
  }
  
  # accounting for age:
  
  if(!is.null(age)){
    ineqs = c();
    pattern = "(>|<)=?"
    for(i in 1:length(age)){
      splt = str_split(age[[i]],pattern=pattern)[[1]]
      ineqsigns = str_extract_all(age[[i]],pattern=pattern)[[1]]
      lsplt = length(splt)
      for(j in 1:(lsplt-1)){
        ineqs = append(ineqs,paste0(splt[j],ineqsigns[j],splt[j+1]))
      }
    
    }
    
    lpm = length(pm)
    pmagemat = matrix(rep(0,lpm^2),nrow=lpm,ncol=lpm)
    pmids = c()
    for(i in pm){
      pmids=append(pmids,i$ID)
    }
    for(i in 1:length(ineqs)){
      splt = str_split(ineqs[[i]],pattern=pattern)
      ineqsign = str_extract(ineqs[[i]],pattern=pattern)
      matches = match(splt[[1]],pmids)
      if(ineqsign=="<"|ineqsign=="<="){
        pmagemat[matches[2],matches[1]] = 1
      }
      else{
        pmagemat[matches[1],matches[2]] = 1
      }
    }
    
    amids = c()
    for(i in am){
      for(j in i$ID){
        amids=append(amids,j)
      }
    }
    lam = length(amids)
    amagemat = matrix(rep(0,lam^2),nrow=lam,ncol=lam)
    for(i in 1:length(am)){
      a = am[[i]]
      lai = length(a$ID)
      minimat = matrix(rep(0,lai^2),lai,lai)
      for(j in 1:lai){
        if(a$FIDX[j]!=0){
          minimat[a$FIDX[j],j] = 1;
        }
        if(a$MIDX[j]!=0){
          minimat[a$MIDX[j],j] = 1;
        }
      }
      inds = match(a$ID,amids)
      amagemat[inds,inds]=minimat
    }
    
    for(i in 1:lpm){
      done = rep(FALSE,lpm)
      changed = TRUE;
      while(changed){
        changed = FALSE;
        for(j in 1:lpm){
          if(pmagemat[i,j]&!done[j]){
            pmagemat[i,] = pmagemat[i,]|pmagemat[j,]
            changed = TRUE;
            done[j] = TRUE;
          }
        }
      }
    }
    for(i in 1:lam){
      done = rep(FALSE,lam)
      changed = TRUE;
      while(changed){
        changed = FALSE;
        for(j in 1:lam){
          if(amagemat[i,j]&!done[j]){
            amagemat[i,] = amagemat[i,]|amagemat[j,]
            changed = TRUE;
            done[j] = TRUE;
          }
        }
      }
    }
    pmagemat = pmagemat-t(pmagemat)
    amagemat = amagemat-t(amagemat)
    
    todo = rep(TRUE,lpm)
    queue = rep(FALSE,lpm)
    pm.order = c();
    am.order = c();
    pmboolean = pmagemat!=0
    amboolean = amagemat!=0
    
    while(sum(todo)){
      if(sum(queue)){
        if(sum(queue)==1){nxt = (1:lpm)[queue]}
        else{
          nxt = max.col(t(apply(pmboolean[,queue],2,sum)));
          nxt = (1:lpm)[todo][nxt]
        }
        pm.order = append(pm.order,nxt);
        todo[nxt] = FALSE;
        queue = (queue|pmboolean[nxt,])&todo
      }
      
      else{
        if(sum(todo)==1) nxt = (1:lpm)[todo]
        else{
          nxt = max.col(t(apply(pmboolean[,todo],2,sum)));
          nxt = (1:lpm)[todo][nxt]
        }
        pm.order = append(pm.order,nxt);
        todo[nxt] = FALSE;
        queue = pmboolean[nxt,]&todo;
      }
    }
    todo = rep(TRUE,lam)
    queue = rep(FALSE,lam)
    while(sum(todo)){
      if(sum(queue)){
        if(sum(queue)==1){nxt = (1:lam)[queue]}
        else{
          nxt = max.col(t(apply(amboolean[,queue],2,sum)));
          nxt = (1:lam)[todo][nxt]
        }
        am.order = append(am.order,nxt);
        todo[nxt] = FALSE;
        queue = (queue|amboolean[nxt,])&todo
      }
      
      else{
        if(sum(todo)==1) nxt = (1:lam)[todo]
        else{
          nxt = max.col(t(apply(amboolean[,todo],2,sum)));
          nxt = (1:lam)[todo][nxt]
        }
        am.order = append(am.order,nxt);
        todo[nxt] = FALSE;
        queue = amboolean[nxt,]&todo;
      }
    }
    pmagemat = pmagemat[pm.order,pm.order]
    amagemat = amagemat[am.order,am.order]
    tmp = matrix(rep(0,(lam+1)^2),nrow=lam+1,ncol=lam+1)
    tmp[1:lam,1:lam]=amagemat
    amagemat = tmp;
    pmr = pm[pm.order]
    amids = amids[am.order]
    pmidz = pmids[pm.order]
    
    pm.index = rep(0,lpm);
    permut = rep(0,lpm);
    position = 1;
    lst.order = lst[pm.order]
    llst = rep(0,lpm)
    for(i in 1:lpm){
      llst[i] = length(lst[[i]])
    }
    full_list = list()
    amids.star = c(amids,"*")
    recurse2 = function(position,permut,pm.index,full_list){
      if(position>1){
        solid = t(amagemat[,permut[1:(position-1)]])*pmagemat[position,1:(position-1)]
        
        if(position>2)solid=apply(solid==-1,2,sum)==0
        else solid = solid!=-1
        solid = solid & c(!((1:lam) %in% permut),TRUE)
      }
      else{
        solid = rep(1,lam+1)
      }
      for(i in 1:llst[position]){
        m = match(lst[[position]][i],amids.star)
        if(solid[m]){
          pm.index[position] = i
          permut[position] = m
          if(position>=lpm){
            new = c();
            for(j in 1:lpm){
              new=c(new,amids.star[permut[j]])
            }
            full_list = append(full_list,list(new))
          }
          else{
            full_list=recurse2(position+1,permut,pm.index,full_list)
          }
        }
      }
      return(full_list)
    }
       
    reslist=recurse2(position,permut,pm.index,full_list)
    res = as.data.frame.matrix(do.call(rbind, reslist))
    names(res)=pmids[pm.order]
    res = res[pmids]
    return(res)
  }
  
  recurse = function(a) {
    if(length(a) == 1)
      return(as.list.default(a[[1]]))
    
    r = lapply(a[[1]], function(val) {
      b = a[-1]
      if(val != "*")
        b = lapply(b, function(vec) vec[vec != val])
      lapply(recurse(b), function(vec) c(val, vec))
    })
    
    unlist(r, recursive = FALSE)
  }
  
  reslist = recurse(lst)
  res = as.data.frame.matrix(do.call(rbind, reslist))
  names(res) = nms
  
  res
}

