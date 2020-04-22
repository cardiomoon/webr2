#' Check for need for transformation
#' @param x A numeric vector
#' @param family	Character The quoted name of a family of transformations.
#' @return A logical value
#' @importFrom car powerTransform
#' @export
#' @examples
#' needTransform(mtcars$mpg)
#' needTransform(iris$Sepal.Length)
needTransform=function(x,family="bcnPower"){
  if(!is.numeric(x)) return(FALSE)
  res=car::powerTransform(x,family=family)
  temp=as.numeric(gsub("< ","",summary(res)$tests$pval[2]))
  temp<0.05
}

#' Autologous transformation of data.frame using bestNormalize
#' @param df A data.frame
#' @param vars2transform Character vector Variable names to transform or NULL
#' @importFrom purrr map_lgl
#' @export
#' @examples
#' autoTransformDf(mtcars)
#' res=autoTransformDf(df=mtcars)
autoTransformDf=function(df,vars2transform=NULL){
  BN=list()
  df2=list()
  if(is.null(vars2transform)) {
     res=map_lgl(df,needTransform)
  } else{
     res= names(df) %in% vars2transform
  }
  for(i in 1:length(res)){
    if(res[i]){
      BN[[names(df)[i]]]=bestNormalize(df[[i]], r = 1, k = 5, warn = F)
      df2[[names(df)[i]]]=seq(min(df[[i]],na.rm=TRUE),max(df[[i]],na.rm=TRUE),length.out=100)
      df[[paste0(names(df)[i],".t")]]=BN[[names(df)[i]]]$x.t
    }else{
      if(is.numeric(df[[i]])){
        df2[[names(df)[i]]]=seq(min(df[[i]],na.rm=TRUE),max(df[[i]],na.rm=TRUE),length.out=100)
      } else{
        df2[[names(df)[i]]]=names(which.max(table(df[[i]])))
      }
      df[[paste0(names(df)[i],".t")]]=df[[names(df)[i]]]
    }
  }
  df2=as.data.frame(df2,stringsAsFactors = FALSE)
  result=list(df=df,
       dfnew=df2,
       BN=BN)
  class(result)="BNdf"
  result
}


#' S3 method of class BNdf
#' @param x An object of class BNdf
#' @param ... Further arguments to be passed to plot_BNdf
#' @export
#' @examples
#' x=autoTransformDf(iris)
#' plot(x)
plot.BNdf=function(x,...){
   plot_BNdf(x,...)
}


#' Plot function for an object of class BNdf
#' @param x An object of class BNdf
#' @param select Choices of BN to plot
#' @importFrom graphics par plot
#' @export
#' @examples
#' x=autoTransformDf(iris)
#' plot_BNdf(x)
plot_BNdf=function(x,select=NULL){
    no=length(x$BN)
    if(is.null(select)) select=1:no
    if(length(select)>1){
      ncol=length(select)%/%2+ifelse(length(select)%%2==1,1,0)
      par(mfrow=c(2,ncol))
    }
    for(i in 1:length(select)){
           plot(x$BN[[select[i]]],main=names(x$BN)[select[i]])
    }
    if(length(select)>1){
      par(mfrow=c(1,1))
    }
}

#' BoxPlot function for an object of class BNdf
#' @param x An object of class BNdf
#' @param select Choices of BN to plot
#' @importFrom ggplot2 facet_grid scale_x_discrete
#' @importFrom tidyselect matches
#' @export
#' @examples
#' x=autoTransformDf(iris)
#' plot_BNdf2(x)
plot_BNdf2=function(x,select=NULL){

  no=length(x$BN)
  if(is.null(select)) select=1:no

  map2_dfr(x$BN[select],names(x$BN[select]),function(x,y){
    df<-x$oos_preds
    df$var=y
    df
  }) %>%
    pivot_longer(cols=!matches("var")) %>%
    ggplot(aes_string(x="name",y="value"))+
    geom_boxplot()+
    facet_grid(var~.,scales="free")+
    labs(x="",y="")+
    theme_bw()+
    scale_x_discrete(guide=guide_axis(n.dodge=2))

}



#' Histogram function for an object of class BNdf
#' @param x An object of class BNdf
#' @param select Choices of BN to plot
#' @importFrom MASS truehist
#' @export
#' @examples
#' x=autoTransformDf(iris)
#' plot_BNdf3(x)
plot_BNdf3=function(x,select=NULL){

  no=length(x$BN)
  if(is.null(select)) select=1:no

  par(mfrow=c(length(select),2))
  for(i in 1:length(select)){
     MASS::truehist(x$BN[[select[i]]]$x,xlab=names(x$BN)[select[i]])
     MASS::truehist(x$BN[[select[i]]]$x.t,xlab=paste0(class(x$BN[[select[i]]]$chosen_transform)[1]," of ",names(x$BN)[select[i]]))
  }
  par(mfrow=c(1,1))

}

#' Predict plot function for an object of class BNdf
#' @param x An object of class BNdf
#' @param select Choices of BN to plot
#' @export
#' @examples
#' x=autoTransformDf(iris)
#' plot_BNdf4(x)
plot_BNdf4=function(x,select=NULL){

  no=length(x$BN)
  if(is.null(select)) select=1:no

  par(mfrow=c(length(select),2))
  for(i in 1:length(select)){
    temp=names(x$BN)[select[i]]
    MASS::truehist(x$BN[[select[i]]]$x.t,xlab=paste0(class(x$BN[[select[i]]]$chosen_transform)[1]," of ",temp))
      plot(x$dfnew[[temp]],predict(x$BN[[select[i]]],newdata=x$dfnew[[temp]]),type="l",col=1,
           main="Best Normalizing Transformation",ylab="g(x)",xlab=temp)
  }
  par(mfrow=c(1,1))

}




#' Fit an autotransformed linear model
#' @param fit An object of class lm
#' @param vars2transform Character vector Variable names to transform or NULL
#' @importFrom stats lm as.formula
#' @importFrom mgcv gam
#' @importFrom purrr pmap_dfc
#' @export
#' @examples
#' data(acs,package="moonBook")
#' fit<- lm(EF ~ weight+TG,data=acs)
#' fit<- lm(EF ~ age+TG,data=acs)  ## age donot need transformation
#' fit<- lm(age ~ EF+TG,data=acs)  ## age donot need transformation
#' fit<- lm(EF ~ age+TG+sex+DM,data=acs)  ## age, sex donot need transformation
#' result=autoTransformFit(fit)
#' fit<- lm(mpg ~ hp+wt,data=mtcars)
#' result=autoTransformFit(fit)
#' plot(result,fill="red")
#' data("autotrader",package="bestNormalize")
#' autotrader$yearsold <- 2017 - autotrader$Year
#' fit<- lm(price ~ mileage + yearsold+status, data = autotrader)
#' result=autoTransformFit(fit)
autoTransformFit=function(fit,vars2transform=NULL){
  # vars2transform=NULL
  df=fit$model
  yvar=names(fit$model)[1]
  xvars=names(fit$model)[-1]
  res=autoTransformDf(df,vars2transform)
  df=res$df

  callstr=deparse(fit$terms)
  gamstr=callstr
  temp=colnames(res$dfnew)
  temp
  for(i in 1:length(temp)){
    if(is.null(vars2transform)){
      if(needTransform(df[[temp[i]]])){
       callstr=gsub(temp[i],paste0(temp[i],".t"),callstr)
      }
    } else{
      if(temp[i] %in% vars2transform){
        callstr=gsub(temp[i],paste0(temp[i],".t"),callstr)
      }
    }
  }

    gamstr=gsub(" ","",gamstr)
    var=unlist(strsplit(gamstr,"~"))[2]
    var=unlist(strsplit(var,"\\+"))
    tevar=var[grepl("\\*",var)]
    tevar
    gsub("*",",",tevar,fixed=TRUE)
    paste0("te(",gsub("*",",",tevar,fixed=TRUE),")")
    xvar=var[!grepl("\\*",var)]
    if(length(tevar)>0){
         gamstr=gsub(tevar,paste0("te(",gsub("*",",",tevar,fixed=TRUE),")"),gamstr,fixed=TRUE)
    }
    if(length(xvar)>0){
      for(i in seq_along(xvar)){
      if(is.numeric(res$dfnew[[xvar[i]]])){
      gamstr=gsub(xvar[i],paste0("s(",xvar[i],")"),gamstr)
      }
    }
  }
  gamstr
  newfit=eval(parse(text=paste0("lm(",callstr,",data=df)")))
  gamfit=gam(as.formula(gamstr),data=df,method="REML")

  newdf=newfit$model
  temp=names(newdf)[-1]

  df3=pmap_dfc(list(a=res$dfnew,b=names(res$dfnew),c=df[names(res$dfnew)]),
              function(a,b,c){
                if(is.null(vars2transform)){
                    if(needTransform(c)){
                        predict(res$BN[[b]],a)
                    } else {
                      a
                    }
                } else{
                    if(b %in% vars2transform){
                      predict(res$BN[[b]],a)
                    } else{
                      a
                    }
                }
  })

  pdata=list()
  p=list()
i=3
  for(i in 1:length(temp)){
    newdata=list()
    gamnewdata=list()
    cateVars=c()
    temp2=gsub("\\.t","",temp[i])
    temp2
    if(is.numeric(df3[[temp2]])){
       newdata[[temp[i]]]=df3[[temp2]]
       gamnewdata[[temp2]]=res$dfnew[[temp2]]
    } else{
       newdata[[temp[i]]]=unique(fit$model[[temp2]])
       gamnewdata[[temp2]]=unique(fit$model[[temp2]])
    }

    for(j in 1:length(temp)){
      if(i==j) next
      temp3=gsub("\\.t","",temp[j])
      if(is.numeric(newdf[[temp[j]]])) {
        newdata[[temp[j]]]=mean(newdf[[temp[j]]],na.rm=TRUE)

        gamnewdata[[temp3]]=mean(res$df[[temp3]],na.rm=TRUE)
      } else{
        # if(temp3==byvar){
        newdata[[temp[j]]]=unique(newdf[[temp[j]]])
        gamnewdata[[temp3]]=unique(res$df[[temp3]])
        cateVars=c(cateVars,temp3)
        # } else{
        # newdata[[temp[j]]]=names(which.max(table(newdf[[temp[j]]])))
        # gamnewdata[[temp3]]=names(which.max(table(res$df[[temp3]])))
        # }
      }
    }

    newdata=as.data.frame(expand.grid(newdata,stringsAsFactors = FALSE),stringsAsFactors = FALSE)
    gamnewdata=as.data.frame(expand.grid(gamnewdata,stringsAsFactors = FALSE),stringsAsFactors = FALSE)
    newdata
    gamnewdata
    yhat=predict(newfit,newdata=newdata,se.fit=TRUE)
    y=yhat$fit
    ymax=yhat$fit+1.96*yhat$se.fit
    ymin=yhat$fit-1.96*yhat$se.fit
    if(is.null(vars2transform)){
      if(needTransform(fit$model[[yvar]])){
        newy=predict(res$BN[[yvar]], newdata = y, inverse = TRUE)
        newymax=predict(res$BN[[yvar]], newdata = ymax, inverse = TRUE)
        newymin=predict(res$BN[[yvar]], newdata = ymin, inverse = TRUE)
      } else{
        newy=y
        newymax=ymax
        newymin=ymin
      }
    } else{
      if(yvar %in% vars2transform){
        newy=predict(res$BN[[yvar]], newdata = y, inverse = TRUE)
        newymax=predict(res$BN[[yvar]], newdata = ymax, inverse = TRUE)
        newymin=predict(res$BN[[yvar]], newdata = ymin, inverse = TRUE)
      } else{
        newy=y
        newymax=ymax
        newymin=ymin
      }
    }

    newdata
    temp[i]
    if(is.numeric(fit$model[[temp2]])){
      df=data.frame(x=res$dfnew[[temp2]],y=newy,ymin=newymin,ymax=newymax,stringsAsFactors = FALSE)
    } else{
      df=data.frame(x=newdata[[temp[i]]],y=newy,ymin=newymin,ymax=newymax,stringsAsFactors = FALSE)
    }
    df$method="Transformed linear fit"
    df
    cateVars
    for(k in seq_along(cateVars)){
        df[[cateVars[k]]]=newdata[[cateVars[k]]]
    }
    predgam <- predict(gamfit, newdata = gamnewdata,se.fit=TRUE)
    pgammax=predgam$fit+1.96*predgam$se.fit
    pgammin=predgam$fit-1.96*predgam$se.fit
    if(is.numeric(fit$model[[temp2]])){
      df2=data.frame(x=res$dfnew[[temp2]],y=predgam$fit,ymin=pgammin,ymax=pgammax)
    } else{
    df2=data.frame(x=gamnewdata[[temp[i]]],y=predgam$fit,ymin=pgammin,ymax=pgammax)
    }
    df2$method="gam fit"

    for(k in seq_along(cateVars)){
      df2[[cateVars[k]]]=gamnewdata[[cateVars[k]]]
    }
    df
    df2

    df=rbind(df,df2)

    pdata[[i]]=df

  }

  result<-list(
      res=res,
      fit=fit,
      newfit=newfit,
      gamfit=gamfit,
      pdata=pdata
  )
  class(result)="autoBN"
  invisible(result)
}



#' S3 method of class autoBN
#' @param x An object of calss autoBN
#' @param ... Further arguments to be passed to plot_autoBN
#' @importFrom ggplot2 geom_pointrange geom_jitter
#' @importFrom cowplot plot_grid
#' @export
#' @examples
#' data(acs,package="moonBook")
#' fit<- lm(EF ~ age+TG+sex,data=acs)  ## age, sex donot need transformation
#' result=autoTransformFit(fit)
#' plot(result,fill="red")
#' plot(result,select=1:2,fill="blue")
#' plot(result,add.gam=TRUE)
#' plot(result,by=sex)
#' plot(result,add.gam=TRUE,by=sex,select=1)
#' plot(result,add.gam=TRUE,by=sex)
#' data("autotrader",package="bestNormalize")
#' autotrader$yearsold <- 2017 - autotrader$Year
#' fit<- lm(price ~ mileage + yearsold+status, data = autotrader)
#' result=autoTransformFit(fit)
#' plot(result,fill="red")
#' plot(result,add.gam=TRUE)
#' plot(result,add.gam=TRUE,select=1)
plot.autoBN=function(x,...){
   plot_autoBN(x,...)
}


#' Plot an object of class autoBN
#' @param x An object of calss autoBN
#' @param add.gam logical Whether or not add gam fit
#' @param show.point logical Whether or not show point
#' @param by Optional name of categorical variable
#' @param select Numeric Plot choices
#' @param ... Further arguments to be passed to geom_ribbon
#' @importFrom ggplot2 geom_pointrange geom_jitter
#' @export
plot_autoBN=function(x,add.gam,show.point,by,select,...){

       # add.gam=TRUE;show.point=TRUE;by=NULL;select=NULL
       # byvar="status"
    byvar=NULL
    colorvar=NULL
    if(missing(add.gam)) add.gam=FALSE
    if(missing(show.point)) show.point=TRUE
    if(!missing(by)) {
      byvar=as.character(substitute(by))
    }
    byvar
    cateVars=names(which(map_lgl(x$fit$model[-1],~(!is.numeric(.x)))))
    cateVars=setdiff(cateVars,byvar)
    cateVars
    count=length(x$pdata)

    if(add.gam) colorvar="method"
    colorvar=c(colorvar,byvar)
    colorvar
    xvars=names(x$fit$model)[-1]
    xvars
    yvar=names(x$fit$model)[1]
    p=list()
    i=3
    for(i in 1:count){
       df=x$pdata[[i]]
       df
       if(!add.gam) df=df[df$method!="gam fit",]

       if(is.numeric(df$x)){
         for(j in seq_along(cateVars)){
           df<-df[df[[cateVars[j]]]==names(which.max(table(x$fit$model[[cateVars[j]]])))[1],]
         }


         p[[i]]<-ggplot(data=df,aes_string(x="x",y="y"))
         if(show.point) p[[i]]<-p[[i]]+geom_point(data=x$fit$model,aes_string(x=xvars[i],y=yvar),alpha=0.1)

         p[[i]]<-p[[i]]+ geom_line(data=df,aes_string(color=colorvar[1]))+labs(x=xvars[i],y=yvar)

         if(length(colorvar)>0) {
           p[[i]]<-p[[i]]+geom_ribbon(data=df,aes_string(fill=colorvar[1],ymax="ymax",ymin="ymin"),alpha=0.6,...)
         } else{
           p[[i]]<-p[[i]]+geom_ribbon(data=df,aes_string(ymax="ymax",ymin="ymin"),alpha=0.6,...)
         }

         if(length(colorvar)>1) p[[i]]<-p[[i]]+facet_wrap(tidyselect::all_of(colorvar[-1]))
         if(length(colorvar)>0) p[[i]]=p[[i]]+theme(legend.position=c(0.8,0.9))


       } else{


          # ggplot(data=df,aes_string(x="x",y="y"))+
          #   geom_pointrange(data=df,aes_string(ymax="ymax",ymin="ymin",color="x"),size=1)

         p[[i]]=ggplot(data=df,aes_string(x="x",y="y"))
         if(show.point){
           p[[i]]=p[[i]]+ geom_jitter(data=x$fit$model,aes_string(x=xvars[i],y=yvar),width=0.2,alpha=0.1)
           # geom_boxplot(data=fit$model,aes_string(x=temp2,y=yvar,fill=colorvar),alpha=0.1,
           #              width=0.2,
           #              position=position_dodge(0.9))
         }
         p[[i]]
         df
         colorvar
         cateVars
         # geom_line(data=df,aes_string(color=colorvar,group=colorvar))+
         # geom_ribbon(data=df,aes_string(fill=colorvar,group=colorvar,ymax="ymax",ymin="ymin"),alpha=0.6)+
         p[[i]]=p[[i]]+geom_pointrange(data=df,aes_string(ymax="ymax",ymin="ymin",color="x"),size=1) +
           labs(x=xvars[i],y=yvar)
         if(length(cateVars)>1) p[[i]]<-p[[i]]+facet_wrap(tidyselect::all_of(cateVars[-1]))
         p[[i]]=p[[i]]+theme(legend.position="none")

       }


    }

    if(!missing(select)) p=p[select]
    if(!is.null(byvar)&(add.gam)) cowplot::plot_grid(plotlist=p)
    else  reduce(p,`+`)
}

