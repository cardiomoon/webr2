# library(PredictABEL)

### reclassification
# data(ExampleData)
# names(ExampleData)[2]="AMD"
#
# form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
# form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
# fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
# fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
# pred1=predict(fit1,type="response")
# pred2=predict(fit2,type="response")
# reclassification(data=ExampleData, cOutcome=2,
#                  predrisk1=pred1, predrisk2=pred2, cutoff <- c(0,.1,.35,1))
#
# ORmultivariate(fit1)
# ORmultivariate(fit2)
# calibrationPlot(y=ExampleData$AMD,pred2)
# boxPlot2(y=ExampleData$AMD,pred2)
# labels=cc("without genetic factors", "with genetic factors")
# predictPlot(list(pred1,pred2),labels=labels)

#' Draw a calibration plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param ... further arguments to be passed to geom_point()
#' @importFrom Hisc::cut2
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' calibrationPlot(y=ExampleData$AMD,pred2)
calibrationPlot=function(y,pred,...){
    y=ExampleData[[2]];pred=pred2
    groups=10
    p<-pred
    matres <- matrix(NA, nrow = groups, ncol = 5)
    sor <- order(p)
    p <- p[sor]
    y <- y[sor]
    groep <- Hmisc::cut2(p, g = groups)
    total <- tapply(y, groep, length)
    predicted <- round(tapply(p, groep, sum), 2)
    observed <- tapply(y, groep, sum)
    meanpred <- round(tapply(p, groep, mean), 3)
    meanobs <- round(tapply(y, groep, mean), 3)
    matres <- cbind(total, meanpred, meanobs, predicted,
                    observed)
    contr <- ((observed - predicted)^2)/(total * meanpred *
                                             (1 - meanpred))
    chisqr <- sum(contr)
    df <- (groups - 2)
    pval <- 1 - pchisq(chisqr, df)
    temp=paste0("atop(chi^2 == ",round(chisqr,2),"(df == ",df,"),italic(p) ")
    if(pval<0.001){
        temp=paste0(temp,"< 0.001)")
    }else {
        temp=paste0(temp,"== ", pval,")")
    }

    df1=as.data.frame(matres)

    ggplot(data=df1,aes(x=meanpred,y=meanobs))+
        geom_segment(x=0,y=0,xend=1,yend=1,color="grey50")+
        geom_point(...)+
        labs(x="Predicted risk",y="Observed risk",title="Calibration plot")+
        xlim(c(0,1))+ylim(c(0,1))+
        annotate("text",x=0.65,y=0.2,label=temp,hjust=0,parse=TRUE)+
        theme_bw()+
        theme(panel.grid = element_blank())
}


#' Draw a box plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param ... further arguments to be passed to geom_boxplot()
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' boxPlot2(y=ExampleData$AMD,pred=pred2)
boxPlot2=function(y,pred,labels,...){
    if(missing(labels)) labels=c("Without disease", "With disease")
    df=data.frame(y=pred2,x=ExampleData$AMD)
    meandiff=mean(df$y[df$x==unique(df$x)[2]])-mean(df$y[df$x==unique(df$x)[1]])
    df$x=factor(df$x,labels=labels)
    ggplot(data=df,aes(x=x,y=y,color=x))+
        geom_boxplot(...)+
        stat_summary(fun="mean",geom="point",pch=15,size=4)+
        labs(x="",y="Predicted risks",title="Difference between groups")+
        annotate("text",x=1.5,y=0.8,
                 label=paste0("Mean diff = ",round(meandiff,3)))+
        theme_bw(base_size = 14)+
        theme(panel.grid = element_blank(),legend.position="None")
}


#' Draw a Predictiveness Plot
#' @param predlist A list of prediction vectors
#' @param labels String Optional labels
#' @param ... Further arguments to be passed to geom_line()
#' @return A ggplot
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' labels=cc("without genetic factors", "with genetic factors")
#' predictPlot(list(pred1,pred2),labels=labels)
predictPlot=function(predlist,labels,size,...){
    count=length(predlist)
    if(missing(labels)) labels <- paste0("Model",1:count)
    if(missing(size)) size=1
    df<-map2_dfr(predlist,1:count,function(pred,i){
        y <- sort(pred)
        n <- length(y)
        x <- (1:n)/n
        df=data.frame(x=x,y=y)
        df$no=i
        df
    })


    df$no=factor(df$no,labels=labels)

    ggplot(data=df,aes(x=x,y=y,color=no,lty=no))+geom_line(size=size,...)+
        labs(x="Cumulative percentage",
             y="Predicted risks",
             title="Predictiveness Curve")+
        theme_bw(base_size = 14)+
        theme(legend.position=c(0.3,0.9),
              legend.title = element_blank(),
              legend.box.background = element_rect(),
              panel.grid = element_blank())
}


#' Draw prior versus posterior plot
#' @param pred1 Numeric vecor as a result of predict.glm
#' @param pred2 Numeric vecor as a result of predict.glm
#' @param xlab,yla,title character
#' @param ... Further arguments to be passed to geom_point()
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' priorPosteriorPlot(pred1,pred2,alpha=0.3,color="red")
priorPosteriorPlot=function(pred1,pred2,xlab,ylab,title,...){
    if(missing(xlab)) xlab="Prior Risk"
    if(missing(ylab)) ylab="Posterior Risk"
    if(missing(title)) title="Prior versus Posterior Risk"
    ggplot(data=NULL)+
        geom_segment(aes(x=0,y=0,xend=1,yend=1),color="black")+
        geom_point(aes(x=pred1,y=pred2),...)+
        labs(x=xlab,y=ylab,title=title)+
        xlim(c(0,1))+ylim(c(0,1))+theme_bw()+
        theme_bw(base_size = 14)+
        theme(panel.grid = element_blank())

}



#' Risk Distribution Plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param labels Character
#' @param ... further arguments to be passed to geom_col()
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' riskDistributionPlot(y=ExampleData$AMD,pred=pred2)
riskDistributionPlot=function(y,pred,labels,...){
    if(missing(labels)) labels=c("Without outcome", "With outcome")

    m <- table(y, cut(pred, seq(0,1.05, 0.05)))
    df=as.data.frame(t(m))
    names(df)=c("x","no","Freq")
    df$x=as.numeric(gsub("\\(|,.*\\]","",df$x))
    df<-df %>% group_by(no) %>% mutate(ratio=Freq/sum(Freq))
    df$no=factor(df$no,labels=labels)
    ggplot(df,aes(x=x,y=ratio,fill=no))+
        geom_col(position="dodge",...)+
        labs(x="Predicted risk",y="Precentage",
             title="Distribution of Predicted Risk")+
        theme_bw()+
        theme(legend.position=c(0.8,0.9),
              legend.title = element_blank(),
              legend.box.background = element_rect(),
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5))
}





#' Draw ROC curves
#' @param data	A data.frame
#' @param cOutcome Column number of the outcome variable.
#' @param predrisk	Vector of predicted risk
#' @param labels Labels given to the ROC curves.
#' @export
#' @examples
#' labels <- c("without genetic factors", "with genetic factors")
#' plotROC2(data=ExampleData, cOutcome=2, predrisk=cbind(pred1,pred2), labels=labels)
plotROC2=function(data, cOutcome, predrisk, labels){
    plotROC(data=data, cOutcome=cOutcome, predrisk=predrisk, labels=labels)
    roc1=roc(data[[cOutcome]],pred1)
    roc2=roc(data[[cOutcome]],pred2)
    auclabel=function(x){
        paste0("AUC: ",round(x$auc,3),"(",round(ci(x)[1],3),"-",round(ci(x)[3],3),")") }
    auclabels=c(auclabel(roc1),auclabel(roc2))

    res=roc.test(roc1,roc2)
    delong="DeLong's test: "
    if(res$p.value<0.001) {
        delong=paste0(delong,"p < 0.001")
    } else {
        delong=paste0(delong,"p =", round(res$p.value,3))
    }
    text(0.7,0.25,label=delong)
    text(0.7,0.20,label=auclabel(roc1))
    text(0.7,0.15,label=auclabel(roc2),col="red")
}




#' Draw ROC curves
#' @param y A numeric vector
#' @param predlist A list of numeric vecors as results of predict.glm
#' @param labels Character
#' @param ... further arguments to be passed to geom_line()
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' labels <- c("without genetic factors", "with genetic factors")
#' ggplotROC(ExampleData$AMD,predlist=list(pred1,pred2),labels=labels)
ggplotROC=function(y,predlist,labels,...) {

if(missing(labels)) labels <- c("Model 1", "Model 2")

count=length(predlist)
temp=list()
df<-predlist %>% map2_dfr(1:count,function(a,b){
    temp[[b]]<<-roc(y~a)
    df=data.frame(y=temp[[b]]$sensitivities,x=1-temp[[b]]$specificities)
    df$no=b
    df
})

auclabel=temp %>% map_chr(function(x){
    paste0("AUC: ",round(x$auc,3),"(",round(ci(x)[1],3),"-",round(ci(x)[3],3),")")
})
df2=data.frame(label=auclabel)
df2$no=1:count
df2$x=0.8
df2$y=0.3-(df2$no/20)
df2$no=factor(df2$no,labels=labels)
res=roc.test(temp[[1]],temp[[2]])

df$no=factor(df$no,labels=labels)
delong="DeLong's test: "
if(res$p.value<0.001) {
    delong=paste0(delong,"p < 0.001")
} else {
    delong=paste0(delong,"p =", round(res$p.value,3))
}
p<-ggplot(data=df)+
    geom_line(aes(x=x,y=y,color=no,group=no,lty=no),size=1,...)+
    labs(x="1-Specificity",y="Sensitivity",
         title="ROC plot")+
    geom_segment(x=0,y=0,xend=1,yend=1,color="grey50",lty=2)+
    theme_bw()+
    theme(legend.position=c(0.8,0.1),
          legend.title = element_blank(),
          legend.box.background = element_rect(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
p+geom_text(data=df2,aes(x=x,y=y,label=label,color=no))+
    annotate("text",x=0.8,y=0.3,label=delong)
}
