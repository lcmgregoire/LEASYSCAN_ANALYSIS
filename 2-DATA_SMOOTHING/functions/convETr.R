convETr <- function(x, y){
  choice <- readline(prompt = "Convert ETr in grams to mm (Y/N): ")
  if(choice=="Y"){
    x[ ,7:ncol(x)] <- y*4/1000
    ETr_F <- x
  }else{
    print("No conversion was done!")
    ETr_F <- x
  } 
  return(ETr_F)
}