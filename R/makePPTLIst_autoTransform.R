#' Make Powerpoint list for autoTransformFit
#' @param fit An object of class lm
#' @param vars2transform Optional names of variables to trnasform
#' @param vanilla logical If true, make a vanilla table
#' @param show.transformDensity,add.gam logical Whether or not add tansformDensity or gam result
#' @export
#' @examples
#' data(autotrader)
#' fit<- lm(price ~ mileage + yearsold +status, data= autotrader )
#' result=makePPTList_autoTransform(fit)
#' result=makePPTList_autoTransform(fit,vars2transform=c('mileage'))
#' result=makePPTList_autoTransform(fit,show.transformDensity=FALSE)
#' result=makePPTList_autoTransform(fit,add.gam=FALSE)
#' result=makePPTList_autoTransform(fit,show.transformDensity=FALSE,add.gam=FALSE)
makePPTList_autoTransform=function(fit,vars2transform=NULL,vanilla=TRUE,show.transformDensity=TRUE,add.gam=TRUE){
     # vars2transform=NULL;vanilla=TRUE;show.transformDensity=TRUE;add.gam=TRUE
    temp=fit2call(fit)
    type="Rcode"
    title="Original Linear model"
    code=paste0("fit=",temp,";fit")
    type=c(type,"Rcode")
    title=c(title,"Summary of Original Linear model")
    code=c(code,"summary(fit)")
    type=c(type,"ggplot")
    title=c(title,"Prediction with Original Linear model")
    code=c(code,"predictModel(fit,alpha=0.05)")

    if(show.transformDensity){
    if(is.null(vars2transform)){
        result=which(unlist(lapply(fit$model,needTransform)))
        result
        for(i in seq_along(result)){
            type=c(type,"ggplot")
            title=c(title,paste0("Best Normalizing Transformation of ",names(result)[i]))
            code=c(code,paste0("transformDensity(fit$model[[",result[i],"]],best=TRUE)"))
        }
    } else{
        for(i in seq_along(vars2transform)){
            type=c(type,"ggplot")
            title=c(title,paste0("Best Normalizing Transformation of ",vars2transform[i]))
            code=c(code,paste0("transformDensity(fit$model[['",vars2transform[i],"']],best=TRUE)"))
        }
    }
    }

    type=c(type,"Rcode")
    title=c(title,"Best Normalizing Result")
    if(is.null(vars2transform)){
        temp=paste0("res<-autoTransformFit(fit);res$res$BN")
        res<-autoTransformFit(fit)
    } else{
        temp=paste0("res<-autoTransformFit(fit,vars2transform=c('",paste0(vars2transform,collapse="','"),"'));res$res$BN")
        res<-autoTransformFit(fit,vars2transform=vars2transform)
    }
    code=c(code,temp)

    type=c(type,"plot","ggplot","plot","plot")
    title=c(title,"Prediction with possible transformation","Boxplot for possible Transformation","Histogram of transformed data",
            "Best Normalizing transformation")
    code=c(code,"plot_BNdf(res$res)","plot_BNdf2(res$res)","plot_BNdf3(res$res)","plot_BNdf4(res$res)")

    type=c(type,"Rcode","ggplot")
    title=c(title,"Summary of Transformed Linear Model","Prediction with Transformed Linear model")
    code=c(code,"summary(res$newfit)","predictModel(res$newfit,alpha=0.05)")

    if(add.gam){
    type=c(type,"Rcode","ggplot")
    title=c(title,"Summary of gam model",
            "Prediction with gam model")
    code=c(code,"summary(res$gamfit)",
           "ggGam::ggGam(res$gamfit)")
    }

    type=c(type,"ggplot")
    title=c(title,"Prediction with Retransformed Data")
    code=c(code,"plot(res,fill='red')")

    for(i in seq_along(res$pdata)){
        type=c(type,"ggplot")
        title=c(title,"Prediction with Retransformed Data")
        temp=paste0("plot(res,add.gam=",add.gam,",select=",i)
        if(!add.gam) temp=paste0(temp,",fill='red'")
        temp=paste0(temp,")")
        code=c(code,temp)
    }

    data.frame(type=type,title=title,code=code,stringsAsFactors = FALSE)
}

