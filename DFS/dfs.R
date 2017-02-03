set.seed(12345)
g = rbinom(16,1,0.2)
graph = matrix(g,4,4)

dfs <- function(node,p,graph){
    nodes_to_visit = setdiff(1:p,node)
    paths = which(graph[,node]>0)
    if(length(paths)==0){
        return(paths)
    } else {
    while(length(nodes_to_visit)>0){
        currentnode = nodes_to_visit[1]
        nodes_to_visit = setdiff(nodes_to_visit,currentnode)
        paths = unique(append(paths,which(graph[,currentnode]>0)))
    }
    return(paths)
}
}

sum(dfs(1,4,graph)==1)
sum(dfs(2,4,graph)==2)
sum(dfs(3,4,graph)==3)
sum(dfs(4,4,graph)==4)
