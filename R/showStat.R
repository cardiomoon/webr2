#' Draw distribution of statistics
#' @param fun distrubution function
#' @param df1,df2 degrees of freedom
#' @param stat numeric The value of statistics
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param alpha numeric confidence level of the interval
#' @param show.area,show.critical,show.point logical
#' @param color character Name of line color
#' @param fill  character Name of fill color
#' @param digits integer indicating the number of decimal places
#' @export
#' @importFrom ggplot2 ggplot stat_function annotate aes geom_hline
#' @examples
#' chisq.test(table(mtcars$am,mtcars$vs))
#' showStat(fun=dchisq,df1=1,stat=0.34754)
#' t.test(mpg~am,data=mtcars)
#' showStat(dt,df1=18.332,stat=-3.7671)
showStat=function(fun=NULL,df1=NULL,df2=NULL,stat=NULL,
                  alternative = c("two.sided", "less", "greater"),
                  alpha=0.05,show.area=TRUE,show.critical=TRUE,
                  show.point=TRUE,color="black",fill="red",digits=3){

    # fun=dt;df1=29;df2=NULL;stat=10;alternative = "two.sided"
    # alpha=0.05;color="black";fill="red";digits=3;show.point=TRUE
    # tempstat="t";show.area=TRUE;show.critical=TRUE
    # fun=dchisq;df1=1;df2=NULL;stat=3.1296
    # alternative="greater"

    alternative=match.arg(alternative)
    if(is.null(fun)) fun=dt
    if(is.null(df1)) df1=5
    tempstat=as.character(substitute(fun))
    tempstat=substr(tempstat,2,nchar(tempstat))
    if(tempstat =="chisq") {alternative="greater"}

    if(alternative=="two.sided") alpha=alpha/2

    p2q=function(p,df1,df2=NULL){
        if(is.null(df2)) temp=paste0("q",tempstat,"(",p,",",df1,")")
        else temp=paste0("q",tempstat,"(",p,",",df1,",",df2,")")
        pos=eval(parse(text=temp))
    }

    if(is.null(df2)) {
        args.list=list(df=df1)
    } else{
        args.list=list(df1=df1,df2=df2)
    }
    xpos1=p2q(1-alpha,df1,df2)
    xpos2=p2q(alpha,df1,df2)
    xmin=ifelse(tempstat=="chisq",0,p2q(0.001,df1,df2))
    xmax=p2q(0.999,df1,df2)

    xmax=max(xmax,stat)
    xmin=min(xmin,stat)
    plotlim=c(xmin,xmax)

    p<-ggplot(NULL,aes(x=plotlim))+
        stat_function(fun=fun,
                      geom="line",
                      color=color,
                      args=args.list)
    if(show.critical) p<-p+stat_function(fun=fun,
                                         geom="area",
                                         fill=fill,
                                         args=args.list,
                                         alpha=0.4,
                                         xlim=c(xpos1,plotlim[2]))


    if((alternative!="greater") & show.critical){
        p<-p+stat_function(fun=fun,
                           geom="area",
                           fill=fill,
                           args=args.list,
                           alpha=0.4,
                           xlim=c(plotlim[1],xpos2))
    }
    p

    stat2q=function(stat,df1,df2=NULL,mode=0){
        if(is.null(df2)) pos=eval(parse(text=paste0("q",tempstat,"(",ifelse(mode,"1-",""),"p",tempstat,"(",stat,",",df1,"),",df1,")")))
        else pos=eval(parse(text=paste0("q",tempstat,"(",ifelse(mode,"1-",""),"p",tempstat,"(",stat,",",df1,",",df2,"),",df1,",",df2,")")))
        pos
    }

    if(!is.null(stat)){
        xpos=stat2q(stat,df1,df2)
        if(xpos>=0) { xlim=c(xpos,xmax)}
        else{ xlim=c(xmin,xpos) }

        if(show.area) p=p+stat_function(fun=fun, geom="area",fill="steelblue",args=args.list,alpha=0.4,xlim=xlim)

        if(alternative!="greater"){
            xpos=stat2q(stat,df1,df2,mode=1)


            if(xpos<0) { xlim=c(xmin,xpos)}
            else{ xlim=c(xpos,xmax) }

            if((abs(stat) < abs(xmin))|(abs(stat) <abs(xmax))) {

                if(show.area) p=p+stat_function(fun=fun, geom="area",fill="steelblue",args=args.list,alpha=0.4,xlim=xlim)
            }
        }
    }

    if(is.null(df2)){
        temp1=paste0("1-p",tempstat,"(",abs(stat),",",df1,")")
    } else{
        temp1=paste0("1-p",tempstat,"(",abs(stat),",",df1,",",df2,")")
    }
    value=eval(parse(text=temp1))
    if(alternative=="two.sided") {
        temp=round(value*2,digits)
    } else{
        temp=round(value,digits)
    }
    if(temp==0){
        temp=paste0("italic(p) < 0.001")
    } else {
        temp=paste0("italic(p) == ",temp)
    }
    hjust=ifelse(stat>0,-0.1,1.1)

    if(abs(xmax-abs(stat))<2) hjust=0.5

    if(is.null(df2)) {
        ypos=fun(stat,df1)
    } else{
        ypos=fun(stat,df1,df2)
    }
    if(ypos<0.01) hjust=ifelse(stat>0,1.1,-0.1)
    p=p+annotate(geom="text",x=stat,y=ypos,label=temp,hjust=hjust,vjust=-0.7,parse=TRUE)+
        annotate(geom="text",x=stat,y=0,label=round(stat,digits),vjust=1.5)
    if(show.point) p=p+annotate(geom="point",x=stat,y=ypos,color="blue")
    p<-p+labs(x=paste0(tempstat," statistic"),y="Probability Density")+theme_bw()+
        theme(panel.grid=element_blank())+
        geom_hline(yintercept=0,color="grey80")
    if(is.null(df2)) label=paste0("df = ",df1)
    else label=paste0("df1 = ",df1,", df2 = ",df2)
    p<-p+annotate("text",x=Inf,y=Inf,label=label,hjust=2,vjust=3)
    p
}



#' Draw distribution of statistic for an object of lm
#' @param fit An object of class lm
#' @importFrom purrr reduce
#' @importFrom stats df dt
#' @export
#' @examples
#' fit=lm(mpg~hp*wt+vs,data=mtcars)
#' showStat_lm(fit)
showStat_lm=function(fit){

    df3=summary(fit)$df[2]
    p=list()
    for(i in 1:nrow(summary(fit)$coef)){
        p[[i]]<-showStat(fun=dt,df1=df3,stat=summary(fit)$coef[i,3])+labs(title=names(fit$coef)[i])
    }
    x=summary(fit)$fstatistic
    p[[i+1]]<-showStat(fun=df,df1=x[2],df2=x[3],stat=x[1],digits=4,alternative="greater")+labs(title=deparse(fit$call))
    purrr::reduce(p,`+`)
}

