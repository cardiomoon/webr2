#' Make new data for predict
#' @param model An object
#' @param length numeric length of continuous variable to to predict
#' @param by character optional factor variable
#' @param type character type argument to be passed to predict.gam
#' @importFrom stats plogis predict
#' @importFrom predict3d restoreData restoreData2 restoreData3
#' @export
#' @examples
#' model=lm(mpg~wt+hp+am,data=mtcars)
#' makeNewData(model)
makeNewData=function(model,length=100,by=NULL,type="response"){

     # length=100;by="am";type="response"
     xvars=names(model$model)[-1]
     if(!is.null(by)) {
          xvars2=setdiff(xvars,by)
     } else{
          xvars2=xvars
     }

     for(j in 1:length(xvars)){
          result=list()

          for(i in 1:length(xvars)){
               x=model$model[[xvars[i]]]

               if(is.factor(x)) {
                    if(i==j) {
                         result[[i]]=levels(x)
                    } else if(xvars[i] %in% by){
                         result[[i]]=levels(x)
                    } else{
                         result[[i]]=levels(x)[1]
                    }
               } else {  ## numeric

                    if(xvars[i] %in% by){
                         result[[i]]=unique(x)
                    } else if(i==j){
                         result[[i]]=seq(from=min(x,na.rm=TRUE),to=max(x,na.rm=TRUE),length.out=length)
                    } else{
                         result[[i]]=mean(x,na.rm=TRUE)
                    }
               }
          }
          result
          names(result)=xvars
          newdata=as.data.frame(expand.grid(result))
          newdata=restoreData(newdata)
          newdata=restoreData2(newdata)
          newdata=restoreData3(newdata)
          names(newdata)


          if ("glm" %in% class(model)) {
            if (model$family$family == "binomial") {
          df1 = as.data.frame(predict(model, newdata = newdata,
                                      type = "link", se.fit = TRUE))
          df1$ymax = df1$fit + 1.96 * df1$se.fit
          df1$ymin = df1$fit - 1.96 * df1$se.fit
          if (type == "response")
            df1[] = lapply(df1, plogis)
            }

          } else{

               df1=as.data.frame(predict(model,newdata=newdata,type=type,se.fit=TRUE))
               df1$ymax=df1$fit+1.96*df1$se.fit
               df1$ymin=df1$fit-1.96*df1$se.fit

          }


          df=cbind(newdata,df1)
          df$xvar=xvars[j]
          if(j==1) {
               final=df
          } else{
               final=rbind(final,df)
          }
     }
     final
}


#' Draw a ggplot with an object of class gam
#' @param model An object of class gam
#' @param select numeric Choices of dependent variables to plot
#' @param point logical Whether or not draw point
#' @param se logical Whether or not draw confidence interval
#' @param by NULL or character optional name of factor variable
#' @param type character type argument to be passed to predict.gam
#' @param interactive logical If true,
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_ribbon ylab facet_wrap xlab scale_x_continuous
#' @importFrom purrr reduce
#' @importFrom ggiraph girafe geom_point_interactive geom_pointrange_interactive opts_hover opts_tooltip
#' @importFrom patchwork patchworkGrob
#' @examples
#' model=lm(mpg~wt+hp+am,data=mtcars)
#' predictModel(model)
#' predictModel(model,interactive=TRUE)
#' predictModel(model,by=am)
#' predictModel(model,by=am,interactive=TRUE)
#' \donttest{
#' model=lm(Sepal.Length~Sepal.Width+Petal.Length+Species,data=iris)
#' summary(model)
#' predictModel(model)
#' predictModel(model,interactive=TRUE)
#' predictModel(model,by=Species)
#' data(colon,package="survival")
#' model=glm(status~rx+age+nodes,data=colon,family=binomial)
#' predictModel(model)
#' predictModel(model,by=rx)
#' predictModel(model,by=rx,type="link")
#' }
predictModel=function(model,select=NULL,point=TRUE,se=TRUE,by,type="response",interactive=FALSE){
     # select=NULL;point=TRUE;se=TRUE;type="response";interactive=FALSE
     byvar=NULL
     if(!missing(by)) byvar=as.character(substitute(by))

     # if(is.null(type)) {
     #      if(model$family$family %in% c("Cox PH","binomial")) {
     #           type="link"
     #      } else{
     #           type="response"
     #      }
     #
     # }
     df1=makeNewData(model,by=byvar,type=type)
    df1
     xvars=names(model$model)[-1]
     xvars2=setdiff(xvars,byvar)
     xvars2
     yvar=names(model$model)[1]
     df=model$model
     df$tooltip=rownames(df)
     p=list()
     fillvar=byvar
     df=restoreData(df)
     df
     if(!is.null(byvar)) {
          df1[[byvar]]=factor(df1[[byvar]])
          df[[byvar]]=factor(df[[byvar]])
     }
     i=1
     for(i in 1:length(xvars2)){
          if(is.null(fillvar)){
               p[[i]]<-ggplot(data=df1[df1$xvar==xvars2[i],])
          }else{
               p[[i]]<-ggplot(data=df1[df1$xvar==xvars2[i],],aes_string(fill=fillvar,group=fillvar))
          }

          if(is.numeric(df[[xvars2[i]]])){

               if(is.null(fillvar)){
                    p[[i]]<-p[[i]]+geom_line(aes_string(x=xvars2[i],y="fit"))
                    if(se) p[[i]]<-p[[i]]+
                         geom_ribbon(aes_string(x=xvars2[i],ymax="ymax",ymin="ymin"),
                                     fill="red",alpha=0.4)
                    if(point) p[[i]]<-p[[i]]+
                    geom_point_interactive(data=df,aes_string(x=xvars2[i],y=yvar,tooltip="tooltip",data_id="tooltip"))
               } else{
                    p[[i]]<-p[[i]]+
                         geom_line(aes_string(x=xvars2[i],y="fit",color=fillvar))

                    if(se) p[[i]]<-p[[i]]+
                         geom_ribbon(aes_string(x=xvars2[i],ymax="ymax",ymin="ymin",fill=fillvar),alpha=0.4)
                    if(point) p[[i]]<-p[[i]]+
                         geom_point_interactive(data=df,aes_string(x=xvars2[i],y=yvar,color=fillvar,tooltip="tooltip",data_id="tooltip"))
               }
          } else{
               p[[i]]<-p[[i]]+
                    geom_pointrange_interactive(aes_string(x=xvars2[i],y="fit",ymax="ymax",ymin="ymin",tooltip=xvars2[i]))
          }


          if(i!=length(xvars2)) {
               p[[i]]=p[[i]]+theme(legend.position="none")
          }
          if(i==1){
               if(("glm" %in% class(model)) & (type=="link")) {
                 p[[i]]<-p[[i]]+ylab("Odds Ratio")
               } else{
               p[[i]]<-p[[i]]+ylab(yvar)
               }
          } else{
               p[[i]]=p[[i]]+ylab("")
          }

     }
     p<-reduce(p,`+`)
     if(interactive){
       girafe(code=print(p),options = list(opts_hover(css = "fill:red;r:3pt;"),opts_tooltip(offx = 10, offy = 10)))
     } else{
       p
     }

}

